# Encoding: utf-8
""" Module for registration of a SpineImaging session
The main steps are:
    1) Cluster the data to get an initial reference volume
    2) Each plane is clustered and xy registered to the initial reference volume from step 1
    3) For each time point, each plane is replaced by the mean xy aligned group it belongs to from step 2
    4) Cluster the data set created in step 4 and estimate the Z position of the resulting groups
    5) Combine the groups from step 4 using the estimated Z position into 4-5 groups
    6) Assign a Z position to each time point in the original data based on the clustering in step 4
    7) xy align the original data to the appropriate reference volume from step 5 based on assignments from step 6

Registration of groups in stage 2 and time points in stage 7:
    1) Compares plane to plane not volume to volume
    2) Use a reference volume the has been expanded to be bigger then the original data
    3) Use knowledge of the possible velocities the brain moves to estimate the likelihood of a given shift
    4) Have priors based on cross-correlation of volumes
    5) Have optimization procedures that try to minimize a metric based on the fano factor of the result
"""

from __future__ import print_function, division

import copy
import logging
import os
import sys
import time
from random import sample

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objs as go
import plotly.offline as py
import plotly.tools as tls
import thunder as td
from pySparkUtils.utils import change, balanced_repartition, fallback
from pySparkUtils.SVD import getSVD
from scipy import signal
from scipy.interpolate import PchipInterpolator
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.morphology import white_tophat
from scipy.optimize import differential_evolution
from scipy.spatial.distance import cdist
from scipy.stats import norm
from scipy.stats.mstats import zscore
from skimage.feature import match_template
from skimage.restoration.inpaint import inpaint_biharmonic
from sklearn.cluster import KMeans
from tsp_solver.greedy_numpy import solve_tsp

from prep.Embedding import getExampleVol, initEmbedDict, prepareCoordinates, getFieldTFormPar
from prep.IO import writeTiff, reloadClean
from prep.Log import add_logging_to_file
from prep.Steps import cropDictStep
from prep.Utils import getTarget, registerByPlane, nanMeanByIndex, getCrop, getBestCrop, log_progress

# Setup logging
logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)

try:
    from itertools import izip
except ImportError:
    # python 3
    izip = zip


def initRegDict(session, data, folderName='intermediates', auPerPhoton=90.0, structureSize=6):
    """ Initialize a dictionary with default values for registration

    :param session: current SpineSession object
    :param data: clean data Images object
    :param folderName: sub-folder name
    :param auPerPhoton: arbitrary value estimating a single photon
    :param structureSize:structure size in pixels in x, y for topHat filtering
    :return: Registration dictionary with default values, creates a directory for intermediate files, starts logging to
    file
    """

    regDict = dict()
    other_vars = locals()
    del other_vars['session']
    regDict.update(other_vars)

    # from session or computed
    regDict['dims'] = data.first().shape
    regDict['path'] = session.path
    regDict['sampleRate'] = float(session.volRate)
    regDict['fullPath'] = os.path.join(session.path + folderName, '')
    regDict['sigmaTime'] = 1000.0 / session.volRate / 2
    regDict['velPerPixXY'] = session.pixSizeXY / (1000 / regDict['sampleRate'])
    regDict['velPerPixZ'] = session.pixSizeZ / (1000 / regDict['sampleRate'])
    regDict['flyLines'] = session.flyLines
    regDict['pixSizeXY'] = session.pixSizeXY
    regDict['pixSizeZ'] = session.pixSizeZ
    regDict['fieldMaskBool'] = session.fieldMaskBool
    regDict['fieldDur'] = float(session.ySize + session.optFlyLines) / session.lineRate * 1000

    if not os.path.isdir(regDict['fullPath']):
        os.mkdir(os.path.dirname(regDict['fullPath']))
    else:
        logger.info('Folder: ' + regDict['fullPath'] + ' exists, files will be overwritten')
    add_logging_to_file('prep.Registration', os.path.join(regDict['fullPath'], 'registration.log'))
    return regDict


@fallback
def change_sc(sc):
    """ prevents spark context from dying when transmitting messages and failing quickly for getFinalGroupsSVD

    :param sc: spark context
    :return: spark context
    """
    return change(sc, spark_rpc_message_maxSize='250', spark_executor_heartbeatInterval='360s',
                  spark_network_timeout='720s', spark_task_maxFailures='1')


def getGlobalSignal(sc, regDict, cropDict=None, minZscore=2.0, minDFF=1.0, regWindow=3000):
    """Estimates global events based on a crop given by cropDict, if None will take all the data.

    :param sc: Spark context
    :param regDict: registration dict
    :param cropDict: cropping dictionary from getCrop function
    :param minZscore: minimal Z score above which to be considered AP (float)
    :param minDFF: minimal Z score above which to be considered AP (float)
    :param regWindow: window in time points to estimate baseline from (int)
    :return:  Adds to regDict:
        1) globalTC: time course of planes that represent global events
        2) putAP: bool array with True for every event above threshold - putative action potential
        3) noAP: Logical not on putAP (all time points without APs)
    """

    data = regDict['data']
    if cropDict is not None:
        data = cropDictStep(data=data, cropDict=cropDict)
    else:
        cropDict = getCrop(data.first(), display=None)
    regDict['cellCrop'] = cropDict

    def getTC(x):
        threshold = np.percentile(x, 95)
        x = x.astype(np.float32)
        x[x < threshold] = np.nan
        return np.nanmean(x)

    globalTC = data.map(getTC).toarray()
    globalTC_BC = sc.broadcast(globalTC)
    regDict['minZscore'] = minZscore
    regDict['minDFF'] = minDFF
    regDict['regWindow'] = regWindow

    def findAP(key):
        """ determine if current time point (key) is an AP

        :param key: current time point
        :return: 1 if AP, 0 otherwise
        """
        from scipy.stats import gaussian_kde
        key = key[0]
        start = int(max((0, key - regWindow / 2)))
        stop = int(min((len(globalTC_BC.value), key + regWindow / 2)))
        y1 = globalTC_BC.value[start:stop]
        y2 = y1[np.logical_not(np.isnan(y1))]
        if np.any(y2):
            if len(y2) > 100:
                kernel = gaussian_kde(y2)
                low = int(np.round(np.percentile(y2, 25)))
                high = int(np.round(np.percentile(y2, 75)))
                step = (high - low) / 100.
                testRange = low + np.array(range(1, 101)) * step
                estMode = testRange[np.argmax(kernel(testRange))]
            else:
                estMode = np.median(y2)
            y3 = y2[(y2 - estMode) < 0] - estMode
            std = np.std(np.hstack((y3, -y3)))
            zScore = (globalTC_BC.value[key] - estMode) / std
            localPutAP = np.logical_and(zScore > minZscore, (globalTC_BC.value[key] / estMode - 1.) > minDFF)
        else:
            localPutAP = 0
        return localPutAP

    index = np.array(range(len(globalTC))).reshape(-1, 1, 1)
    putAP = td.images.fromarray(index, engine=sc).map(findAP).toarray()
    globalTC_BC.unpersist()
    regDict['globalTC'] = globalTC
    regDict['putAP'] = putAP
    noAP = np.ones(data.shape[0], dtype=bool)
    noAP[putAP] = False
    regDict['noAP'] = noAP
    APNum = np.where(regDict['putAP'])[0].shape[0]
    length = regDict['data'].shape[0]
    meanFR = APNum / (length / regDict['sampleRate'])
    regDict['meanFR'] = meanFR
    # plot
    data = list()
    data.append(dict(y=globalTC + 1000, type='scatter', name='TC'))
    data.append(dict(y=putAP * 1000, type='scatter', name='AP'))
    layout = dict(title='Mean firing rate: %.2fHz' % meanFR)
    py.iplot(dict(data=data, layout=layout))


def getInitClustersSVD(regDict, initClusterNum=30, initIterKMeans=35, initNInit=50, initTol=1e-8, initK=40):
    """ Clusters the data set into initClusterNum clusters by first extracting the first initK features of the data
     using SVD and then clustering using K-Means in feature space

    :param regDict: registration dict
    :param initClusterNum: number of clusters for initial reference selection (int)
    :param initIterKMeans: number of KMeans iterations (int)
    :param initNInit: number of initializations for the K-Means
    :param initTol: Tolerance for the K-Means
    :param initK: number of features to get from the SVD
    :return: adds to regDict:
        1) initClusters: cluster assignment per time point
        2) initClusterNums: number of time points per cluster
    """
    local = locals()
    del local['regDict']
    regDict.update(local)

    # SVD on the data
    t = time.time()
    data = regDict['data']
    indexNoAP = np.where(regDict['noAP'])[0]
    U, _, _ = getSVD(data, initK, getComponents=False, getS=False, normalization='mean')
    current = time.time() - t
    m, s = divmod(current, 60)
    logger.info('SVD:: %02d:%02d' % (m, s))

    # K-Means cluster
    model = KMeans(n_clusters=initClusterNum, init='k-means++', n_init=initNInit,
                   max_iter=initIterKMeans, tol=initTol, precompute_distances='auto',
                   verbose=0, random_state=1, copy_x=False, n_jobs=-1,
                   algorithm='auto').fit(U[indexNoAP, :])
    current = time.time() - t
    m, s = divmod(current, 60)
    logger.info('Cluster:: %02d:%02d' % (m, s))
    clusters = model.predict(U)
    centers = model.cluster_centers_

    # get unique cluster
    uIndex, loc = np.unique(clusters, return_inverse=True)
    clusters2 = np.zeros(clusters.shape)
    for i in range(len(uIndex)):
        clusters2[loc == i] = i
    clusters = clusters2
    centers = centers[uIndex]
    dist = cdist(centers, centers)

    # plot
    ax1 = plt.subplot(1, 1, 1)
    y, x, f = ax1.hist(clusters, bins=len(uIndex), color='b')
    ax1.set_ylabel('Time points', color='b')
    ax1.set_xlabel('Cluster#')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    binCenters = 0.5 * (x[1:] + x[:-1])
    ax2 = ax1.twinx()
    ax2.plot(binCenters, np.sum(dist, axis=0), '*r')
    ax2.set_ylabel('Sum distances', color='r')
    plt.grid(False)
    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    # return data
    regDict['initClusters'] = clusters
    regDict['initClusterNums'] = y


def runInitTarget(sc, regDict, rankList=(1, 2, 3), session=None, TForm=None):
    """ A helper function to take the (rankList) clusters and find a target using cross-correlation to the mean

    :param sc: Spark Context
    :param regDict: registration dict
    :param rankList: iterable of ranks for the clusters to select
    :param session: SpineSession object
    :param TForm: transformation dict to embed targets for manual selection
    :return: TForm
    """
    if TForm is None:
        if session is None:
            raise ValueError('If TForm not provided please provide a session object')
        embedDict = initEmbedDict(session)
        prepareCoordinates(embedDict)
        xStart = session.pipelines['Clean'].steps['crop2'].args['xStart']
        xStop = session.pipelines['Clean'].steps['crop2'].args['xStop']
        logger.info('Creating a TForm with xStart: %d, and xStop %d' % (xStart, xStop))
        TForm = getFieldTFormPar(sc, embedDict, xDim=(xStart, xStop))
        embedDict['TForm1'] = TForm
        session.embedDict = embedDict
    for rank in rankList:
        getInitTarget(sc=sc, regDict=regDict, rank=rank, initCC=90)
        regTargetEmbedded = getExampleVol(sc, data=regDict['regTarget'], TFormDict=TForm, project=True)
        writeTiff(regDict['fullPath'], regTargetEmbedded, 'regTarget' + str(rank))
    return TForm


def getInitTarget(sc, regDict, rank=1, initCC=90):
    """ gets the initial reference volume using the 'rank' biggest cluster followed by taking the initCC most
    correlated points and registering them to the mean

    :param sc: Spark context
    :param regDict: registration dict
    :param rank: which cluster to take 1 - most number of timepoints (int)
    :param initCC: CC cutoff for initial target creation (int)
    :return: adds to regDict:
        1) initTargetIndex: which time points are in the rank highest cluster
        2) initTarget: initial target before registration
        3) regTarget: initial target after registration
        4) initPoints: point used to create initTarget
        5) targetShifts: shifts used to create regTarget
    """
    # get params
    data = regDict['data']
    regDict['initCC'] = initCC
    length = regDict['data'].shape[0]
    y = regDict['initClusterNums']
    clusters = regDict['initClusters']
    maxCluster = np.argpartition(y, -rank)[-rank:][0]
    maxIndex = np.where(clusters == maxCluster)[0]

    # crop out cell body (the crop used to get global signal)
    noCellIndex = np.where(np.logical_not(regDict['cellCrop']['zIndex']))[0]
    if noCellIndex.shape[0] > 0:
        data2 = data[:, :, :, noCellIndex]
    else:
        data2 = data

    # plot where the time points of the selected group came from
    plt.figure()
    plt.plot(maxIndex, np.ones(maxIndex.shape), 'r*')
    plt.xlim(0, length)
    plt.title('Time points taken for initTarget')

    # plot the CC
    initPoints = getTarget(data2[maxIndex, :, :, :], initCC, 1, mode='index')
    initTarget = data[initPoints, :, :, :].mean().toarray()

    # get the reference volume
    targetShifts, regTarget = registerByPlane(sc, data=data[initPoints, :, :, :], target=initTarget, upSample2=10,
                                              doShift=True)
    regTarget = regTarget.mean().toarray()

    # save to regDict
    regDict['initTargetIndex'] = maxIndex
    regDict['initTarget'] = initTarget
    regDict['initPoints'] = initPoints
    regDict['regTarget'] = regTarget
    regDict['targetShifts'] = targetShifts


@fallback
def getAllGroupsSVD(sc, session, regDict, planeClusterNum=800, planeIterKMeans=20, planeNInit=10, planeTol=1e-5,
                    planeK=40, planeIterSVD=50, full_path=None, change_sc_cluster=True):
    """ use SVD to get planeK components and cluster each plane to planeClusterNum that should represent similar
     timepoints in x, y z shifts or activity

    :param sc: Spark context
    :param session: Session object
    :param regDict: registration dictionary
    :param planeClusterNum:  number of clusters per plane
    :param planeIterKMeans: number of iteration for the KMeans clustering
    :param planeNInit: number of random initializations for the KMeans clustering
    :param planeTol: the tolerance for the KMeans clustering
    :param planeK: the number of SVD components
    :param planeIterSVD: The number od SVD iterations
    :return: sc,  adds to regDict:
        1) groups: cluster assignment per time point per plane
        2) allMeanGrpVol: list per plane of a 4d numpy array (groupID, x, y, z) mean volumes
    """
    # todo: make planeClusterNum dependent on length / variability in the data
    # init params
    t = time.time()
    data = regDict['data']
    regDict['planeClusterNum'] = planeClusterNum
    regDict['planeIterKMeans'] = planeIterKMeans
    regDict['planeNInit'] = planeNInit
    regDict['planeTol'] = planeTol
    regDict['planeK'] = planeK
    regDict['planeIterSVD'] = planeIterSVD
    planeNum = data.shape[3]
    # get all SVDs
    if change_sc_cluster:
        sc = change(sc, spark_executor_cores='5', spark_task_cpus='5', spark_executor_memory='30g')
        data = reloadClean(sc, session, full_path=full_path)

    def merge_arrays(kv):
        plane, arrays = kv
        list_arrays = list(arrays)
        sorted_list = sorted(list_arrays, key=lambda x: (x[0][0], x[0][1]))
        sorted_array = np.array(list(map(lambda x: x[1], sorted_list)), dtype=np.int16).T
        return plane, sorted_array

    def bring_key_out(kv):
        key, value = kv
        return key[2], (key[:2], value)

    clean_plane = data.toseries().tordd().map(bring_key_out).groupByKey().map(merge_arrays)

    def get_clusters(kv):
        plane, data = kv
        import os
        os.environ['MKL_NUM_THREADS'] = '5'
        os.environ['OMP_NUM_THREADS'] = '5'
        from sklearn.decomposition import TruncatedSVD
        tsvd = TruncatedSVD(planeK, algorithm="randomized", n_iter=10)
        data_reduced = tsvd.fit_transform(data)
        clusters = KMeans(n_clusters=planeClusterNum, init='k-means++', n_init=planeNInit,
                          max_iter=planeIterKMeans, tol=planeTol, precompute_distances='auto',
                          verbose=0, random_state=1, copy_x=False, n_jobs=-1,
                          algorithm='auto').fit_predict(data_reduced)
        return plane, clusters

    clusters = clean_plane.map(get_clusters).collect()
    sorted_list = sorted(clusters, key=lambda x: x[0])
    groups = np.array(list(map(lambda x: x[1], sorted_list)), dtype=np.int16).T
    if change_sc_cluster:
        sc = change_sc(sc)
        # reload clean
        data = reloadClean(sc, session, name=None, returnRegDict=False, returnTForm=False, full_path=full_path)
        regDict['data'] = data
    allMeanGrpVol = [None] * planeNum
    for i in log_progress(range(planeNum)):
        allMeanGrpVol[i] = nanMeanByIndex(data, groups[:, i], largeGroupMean=False)

    regDict['groups'] = groups
    regDict['allMeanGrpVol'] = allMeanGrpVol
    current = time.time() - t
    m, s = divmod(current, 60)
    logger.info('Done: %02d:%02d' % (m, s))
    return sc


def getExpandedTarget(sc, regDict, nIter, flyLinesFlag=True, expansionFactor=2.0, nanThreshold=200, cutoffCC=50,
                      structureSize=3, inPaint=1, zscoreFlag=False, tophatFlag=True):
    """ Iteratively register the groups from getAllGroups to regTarget in an expanded space in X, Y then crop optimally
    using nanThreshold as maximum number of NaNs in the volume

    :param sc: Spark context
    :param regDict: registration dict
    :param nIter: number of iteration to align and crop the reference volume
    :param flyLinesFlag: crop fly lines after shifting (True)
    :param expansionFactor: factor to expand in x and y (2.0)
    :param nanThreshold: number of continuous NaNs permitted in the expanded volume
    :param cutoffCC: cutoff for CC while creating target
    :param structureSize: size for top hat filter
    :param inPaint: 0 = leave NaNs, 1 = use inPaint, 2 = use nan_to_num
    :param zscoreFlag: flag to do zscore
    :param tophatFlag: flag to do top hat filtering
    :return: sc &
        Adds to regDict:
        1) expandedTarget: the expanded reference volume
        2) grpVolShifts: shifts per group used to create the volume
        3) grpVolShiftsTime: shifts per time point to create the volume
        4) intermediate results in *Dict keys
    """
    # init params
    # todo: nanThreshold should be % pixels per plane
    start = time.time()
    regDict['structureSizeExtendedTarget'] = structureSize
    allMeanGrpVol = regDict['allMeanGrpVol']
    zNum = len(allMeanGrpVol)
    regDict['expansionFactor'] = expansionFactor
    regTarget = copy.deepcopy(regDict['regTarget'])
    dims = regDict['dims']
    regDict['nanThresholdExpandedTarget'] = nanThreshold
    regDict['ExpandedTargetCC'] = cutoffCC
    flyLines = regDict['flyLines']
    refTarget = regTarget
    xMeanLast = 0
    yMeanLast = 0
    reference3DDict = dict()
    expendedDict = dict()
    shiftsDict = dict()
    returnedDict = dict()
    if not isinstance(cutoffCC, list):
        cutoffCC = [cutoffCC]
    while len(cutoffCC) <= nIter:
        cutoffCC.append(cutoffCC[-1])

    # paralleling the volumes per plane per group 5d array
    if 'allMeanGrpVolPar' not in regDict.keys() or 'allMeanGrpVolPar2' not in regDict.keys():
        logger.info('Parallelizing allMeanGrpVolPar')
        allMeanGrpVol2 = np.stack(allMeanGrpVol)
        allMeanGrpVolPar = td.images.fromarray(allMeanGrpVol2.transpose((1, 0, 2, 3, 4)),
                                               npartitions=sc.defaultParallelism * 4, engine=sc)
        allMeanGrpVolPar.cache()
        allMeanGrpVolPar.count()
        regDict['allMeanGrpVolPar'] = allMeanGrpVolPar
        allMeanGrpVolPar2 = td.images.fromarray(allMeanGrpVol2, npartitions=sc.defaultParallelism * 4, engine=sc)
        allMeanGrpVolPar2.cache()
        allMeanGrpVolPar2.count()
        regDict['allMeanGrpVolPar2'] = allMeanGrpVolPar2
    else:
        logger.info('Using existing allMeanGrpVolPar')
        allMeanGrpVolPar = regDict['allMeanGrpVolPar']
        allMeanGrpVolPar2 = regDict['allMeanGrpVolPar2']

    # get shift function
    def shift3D(kv):
        key, array = kv
        shift = np.zeros((zNum, 4))
        ref = copy.deepcopy(refTargetBC.value)
        source2 = copy.deepcopy(array).astype(np.float64)
        for z2 in range(zNum):
            source = copy.deepcopy(source2[z2, :, :, :])
            for j in range(source.shape[2]):
                source[:flyLines[j], :, j] = np.nan
                # source[:, :, j] = inpaint_biharmonic(source[:, :, j], np.isnan(source[:, :, j]))
                source[:, :, j] = np.nan_to_num(source[:, :, j])
                if tophatFlag:
                    source[:, :, j] = white_tophat(source[:, :, j], structureSize)
            if zscoreFlag:
                source = zscore(source, axis=2)
            c = match_template(ref, source, pad_input=True)
            if (maxShift[0] * 2) < c.shape[0]:
                crop = np.ceil(float(c.shape[0]) / 2 - maxShift[0]).astype(int)
                c = c[crop:-crop, :, :]
            if (maxShift[1] * 2) < c.shape[1]:
                crop = np.ceil(float(c.shape[1]) / 2 - maxShift[1]).astype(int)
                c = c[:, crop:-crop, :]
            zCenter = int(np.floor(zNum / 2))
            if maxShift[2] == 0:
                c = c[:, :, zCenter]
            else:
                c = c[:, :, (zCenter - maxShift[2]):(zCenter + maxShift[2] + 1)]
            bestShift = np.unravel_index(np.argmax(c), c.shape)
            bestShift -= np.floor(np.asarray(c.shape) / 2)
            if maxShift[2] == 0:
                bestShift = np.hstack((bestShift, 0.))
            shift[z2, :3] = bestShift
            shift[z2, 3] = np.max(c)
        return shift

    # apply shift function
    def itkRealign3D(kv):
        from SimpleITK import ResampleImageFilter, GetImageFromArray, TranslationTransform, GetArrayFromImage
        key, array = kv
        resampler = ResampleImageFilter()
        CC = shiftsOutBC.value[:, key[0], 3] * -1
        cutoff = np.percentile(CC, cutoffCC[k])
        indexes = np.where(CC > cutoff)[0]
        realigned = np.zeros((indexes.shape[0], int(dims[0] * expansionFactor), int(dims[1] * expansionFactor)))
        for i, groupID in enumerate(indexes):
            current = array[groupID, :, :, key[0]].astype('float32')
            if flyLinesFlag:
                current[:flyLines[key], :] = np.nan
            moving = GetImageFromArray(current)
            transform = TranslationTransform(2, shiftsOutBC.value[groupID, key[0], 1:3])
            resampler.SetTransform(transform)
            resampler.SetSize((int(dims[1] * expansionFactor), int(dims[0] * expansionFactor)))
            resampler.SetDefaultPixelValue(np.NAN)
            resampler.SetOutputOrigin(outOrigin)
            realigned[i, :, :] = GetArrayFromImage(resampler.Execute(moving))
        return np.nanmean(realigned, axis=0)

    for k in range(nIter):
        # get the shifts
        for z in range(refTarget.shape[2]):
            refTarget[:flyLines[z], :, z] = np.nan
            refTarget[:, :, z] = inpaint_biharmonic(refTarget[:, :, z], np.isnan(refTarget[:, :, z]))
            if tophatFlag:
                refTarget[:, :, z] = white_tophat(refTarget[:, :, z], structureSize)
        if zscoreFlag:
            refTarget = zscore(refTarget, axis=2)
        maxShift = np.array([refTarget.shape[0] / 2, refTarget.shape[1] / 2, 0])
        refTargetBC = sc.broadcast(refTarget)
        shifts3D = allMeanGrpVolPar.map(shift3D, with_keys=True).toarray()
        refTargetBC.unpersist()
        shiftsDict[k] = copy.deepcopy(shifts3D)
        # center the shifts
        if k == nIter - 1:
            break
        xMean = (np.percentile(shifts3D[:, :, 0], 98) + np.percentile(shifts3D[:, :, 0], 2)) / 2. - xMeanLast
        yMean = (np.percentile(shifts3D[:, :, 1], 98) + np.percentile(shifts3D[:, :, 1], 2)) / 2. - yMeanLast
        shifts3D[:, :, 0] = shifts3D[:, :, 0] - xMean
        shifts3D[:, :, 1] = shifts3D[:, :, 1] - yMean
        xMeanLast = copy.deepcopy(xMean)
        yMeanLast = copy.deepcopy(yMean)
        # get the new expanded target
        shifts3D[:, :, 2] = 0
        shiftsOut = copy.deepcopy(shifts3D)
        shiftsOut = -shiftsOut[:, :, [2, 1, 0, 3]]
        shiftsOutBC = sc.broadcast(shiftsOut)
        outOrigin = (-dims[1] / 2., -dims[0] / 2.)
        logger.info('Iteration  %d, target shape: %s' % (k + 1, refTarget.shape))
        reference3D = allMeanGrpVolPar2.map(itkRealign3D, with_keys=True).toarray().transpose(1, 2, 0)
        shiftsOutBC.unpersist()
        reference3DDict[k] = reference3D
        expandedTarget = copy.deepcopy(reference3D)
        # iteratively find best target size
        maxShifts = np.abs(shifts3D[:, :, :2]).max(axis=(0, 1))
        returned, expandedTarget = getBestCrop(sc, target=expandedTarget, maxShifts=maxShifts,
                                               nanThreshold=nanThreshold, inPaint=inPaint)
        returnedDict[k] = copy.deepcopy(returned)
        expendedDict[k] = copy.deepcopy(expandedTarget)
        refTarget = copy.deepcopy(expandedTarget)
        logger.info('expended shape: %s' % (refTarget.shape,))
        m, s = divmod(time.time() - start, 60)
        logger.info('Time:  %02d:%02d' % (m, s))

    regDict['expandedTarget'] = expandedTarget
    regDict['grpVolShifts'] = shifts3D[:, :, :3]
    regDict['reference3DDict'] = reference3DDict
    regDict['expandedDict'] = expendedDict
    regDict['shiftsDict'] = shiftsDict
    regDict['returnedDict'] = returnedDict

    # assign the shifts to time points
    groups = (regDict['groups']).astype(int)
    grpVolShiftsTime = np.zeros((groups.shape[0], groups.shape[1], 2), dtype=float)
    for timePoint, allPlaneGroups in enumerate(groups):
        for plane, group in enumerate(allPlaneGroups):
            grpVolShiftsTime[timePoint, plane, :] = -shifts3D[group, plane, :2]
    regDict['grpVolShiftsTime'] = grpVolShiftsTime


def runGroupAlignment(sc, regDict, TForm1, percentile=30):
    """ A helper function to run and output to disk the results of registering the groups per plane using the volume
     shifts or after optimization of per plane registration.

    :param sc: Spark Context
    :param regDict: registration dict
    :param TForm1: TFrom to embed to 3d volume for manual inspection
    :param percentile: for optimization see function XYAlignedGroupsOptimization
    :return: writes to disk the mean volumes after registration
    """
    getXYAlignedGroups(sc, regDict, optimization=False, useVolumeShifts=True)
    allGroups = np.array(regDict['planeGrpAlign'])
    sz = allGroups.shape
    dims = regDict['dims']
    allGroups = allGroups.reshape(sz[0], sz[1], dims[0], dims[1]).transpose(0, 1, 3, 2)
    allGroupsEmbedded = getExampleVol(sc, allGroups.transpose(1, 3, 2, 0), TForm1, project=True)
    writeTiff(regDict['fullPath'], allGroupsEmbedded.transpose(1, 2, 0), 'allGroupsVolume')
    writeTiff(regDict['fullPath'], np.expand_dims(np.nanmean(allGroupsEmbedded, axis=0), axis=2),
              'allGroupsVolumeMean')

    # get % nans and fano factor of the data using the prior to use for optimization as the % change and not
    # absolute values
    regDict['nans_group_volume'] = np.sum(np.isnan(allGroups)).astype(float) / np.prod(sz) * 100
    meanAllGroups = np.nanmean(allGroups, axis=1)
    t = np.nanpercentile(meanAllGroups, 50)
    mask = allGroups < t
    allGroups[mask] = np.nan
    mean = np.nanmean(allGroups, axis=1)
    var = np.nanvar(allGroups, axis=1)
    regDict['fano_group_volume'] = np.nanpercentile(var / mean, percentile)
    logger.info('Group volume fano: %f.2, Nans: %f.2' % (regDict['fano_group_volume'], regDict['nans_group_volume']))

    # optimization
    XYAlignedGroupsOptimization(sc, regDict, percentile=percentile)
    allGroups = np.array(regDict['planeGrpAlign'])
    sz = allGroups.shape
    dims = regDict['dims']
    allGroups = allGroups.reshape(sz[0], sz[1], dims[0], dims[1]).transpose(0, 1, 3, 2)
    allGroupsEmbedded = getExampleVol(sc, allGroups.transpose(1, 3, 2, 0), TForm1, project=True)
    writeTiff(regDict['fullPath'], allGroupsEmbedded.transpose(1, 2, 0), 'allGroupsOpto')
    writeTiff(regDict['fullPath'], np.expand_dims(np.nanmean(allGroupsEmbedded, axis=0), axis=2),
              'allGroupsOptoMean')


def XYAlignedGroupsOptimization(sc, regDict, percentile=30):
    """ optimization of group alignment in x y by minimizing fano factor and penalizing the number of NaNs

    :param sc: Spark context
    :param regDict: registration dict
    :param percentile:
    :return: adds to regDict optimized params:
        1) lambdaMotion2
        2) sigmaTime
        3) sigmaPix
    """
    # init params
    start = time.time()
    fieldMaskBool = regDict['fieldMaskBool']
    fieldDur = regDict['fieldDur']
    fieldTimes = np.asarray(range(0, len(fieldMaskBool))) * fieldDur
    fieldTimes = fieldTimes[fieldMaskBool]
    regDict['fieldTimes'] = fieldTimes
    # creating indexing
    allMeanGrpVol = regDict['allMeanGrpVol']
    si = regDict['expandedTarget'].shape
    totalGroups = sum([x.shape[0] for x in allMeanGrpVol])
    regDict['totalGroups'] = totalGroups

    counter = 0
    parIdx = [None] * totalGroups
    for i in range(0, si[2]):
        for j in range(0, allMeanGrpVol[i].shape[0]):
            parIdx[counter] = (i, j)
            counter += 1
    regDict['parIdx'] = parIdx
    if 'allMeanGrpVolRDD' not in regDict or regDict['allMeanGrpVolRDD'].context._jvm is None:
        logger.info('XYAlignedGroupsOptimization making allMeanGrpVolAll')
        allMeanGrpVolAll = np.concatenate(regDict['allMeanGrpVol'])
        allMeanGrpVolRDD = sc.parallelize(izip(iter(parIdx), iter(allMeanGrpVolAll)),
                                          min(len(parIdx), sc.defaultParallelism * 10))
        allMeanGrpVolRDD.cache()
        allMeanGrpVolRDD.count()
        regDict['allMeanGrpVolRDD'] = allMeanGrpVolRDD

    result = differential_evolution(optimizeFanoGroups, bounds=[(0.05, 50), (5, 25), (10, 800), (5, 6)],
                                    args=(sc, regDict, percentile, start),
                                    maxiter=5, tol=1e-2, popsize=4, mutation=(0.5, 1), disp=True, polish=False)
    regDict['optimizationResultGroups'] = result
    lambdaMotion = result.x[0] / 100.
    sigmaTime = result.x[1]
    sigmaPix = result.x[2] / 100.
    structureSize = int(result.x[3])
    logger.info('Running on all data')
    getXYAlignedGroups(sc, regDict, lambdaMotion=lambdaMotion, sigmaPix=sigmaPix, sigmaTime=sigmaTime,
                       structureSize=structureSize)


def optimizeFanoGroups(params, sc, regDict, percentile, start):
    """ test params with getWeightedShifts by calculating the % change in # of NaNs and fano factor

    :param params: tuple of (lambda_motion, sigma_time, sigma_pix)
    :param sc: Spark Context
    :param regDict: registration dict
    :param percentile: to use for fano factor
    :param start: start time of the optimization function
    :return:
    """
    lambda_motion = params[0] / 100.
    sigma_time = params[1]
    sigma_pix = params[2] / 100.
    structureSizeTemp = int(params[3])
    getXYAlignedGroups(sc, regDict, optimization=True, lambdaMotion=lambda_motion, sigmaTime=sigma_time,
                       sigmaPix=sigma_pix, structureSize=structureSizeTemp)
    allGroups = np.array(regDict['planeGrpAlign'])
    sz = allGroups.shape
    dims = regDict['dims']
    allGroups = allGroups.reshape(sz[0], sz[1], dims[0], dims[1])
    meanAllGroups = np.nanmean(allGroups, axis=1)
    t = np.nanpercentile(meanAllGroups, 50)
    mask = allGroups < t
    sz2 = allGroups.shape
    Nans = np.sum(np.isnan(allGroups)).astype(float) / np.prod(sz2) * 100
    allGroups[mask] = np.nan
    mean = np.nanmean(allGroups, axis=1)
    var = np.nanvar(allGroups, axis=1)
    fano = np.nanpercentile(var / mean, percentile)
    assert fano > 0
    fano_rel = fano / regDict['fano_group_volume'] * 100 - 100
    Nans_rel = Nans - regDict['nans_group_volume']
    resultTemp = fano_rel + Nans_rel
    m, s = divmod(time.time() - start, 60)
    logger.info(
        'L:% 6.3f,T:% 4.1f,P:% 4.2f,S:%d,N:%.2f,f:%.2f,r:%.3f, %02d:%02d' %
        (lambda_motion, sigma_time, sigma_pix, structureSizeTemp, Nans_rel, fano_rel, resultTemp, m, s))
    return resultTemp


def getXYAlignedGroups(sc, regDict, optimization=False, useVolumeShifts=False, lambdaMotion=0.15, microMotionGain=1.0,
                       sigmaPix=1.0, sigmaTime=36, structureSize=3):
    """ register in X and Z the groups

    :param sc: SparkContext
    :param regDict: registration dict
    :param optimization: flag to indicate if called from the optimization function
    :param useVolumeShifts: flag to just shift the data using the priors
    :param lambdaMotion:
    :param microMotionGain:
    :param sigmaPix:
    :param sigmaTime:
    :param structureSize:
    :return: adds to regDict:
        1) planeGrpAlign
        2) shiftsGroups
        3) shiftsGroupsTime
    """
    start = time.time()
    if not optimization:
        fieldMaskBool = regDict['fieldMaskBool']
        fieldDur = regDict['fieldDur']
        fieldTimes = np.asarray(range(0, len(fieldMaskBool))) * fieldDur
        fieldTimes = fieldTimes[fieldMaskBool]
        regDict['fieldTimes'] = fieldTimes
        # creating indexing
        allMeanGrpVol = regDict['allMeanGrpVol']
        si = regDict['expandedTarget'].shape

        totalGroups = sum([x.shape[0] for x in allMeanGrpVol])
        regDict['totalGroups'] = totalGroups

        counter = 0
        parIdx = [None] * totalGroups
        for i in range(0, si[2]):
            for j in range(0, allMeanGrpVol[i].shape[0]):
                parIdx[counter] = (i, j)
                counter += 1
        regDict['parIdx'] = parIdx
        if 'allMeanGrpVolRDD' not in regDict or regDict['allMeanGrpVolRDD'] is None:
            logger.info('getXYAlignedGroups making allMeanGrpVolAll')
            allMeanGrpVolAll = np.concatenate(regDict['allMeanGrpVol'])
            allMeanGrpVolRDD = sc.parallelize(izip(iter(parIdx), iter(allMeanGrpVolAll)),
                                              min(len(parIdx), sc.defaultParallelism * 10))
            allMeanGrpVolRDD.cache()
            allMeanGrpVolRDD.count()
            regDict['allMeanGrpVolRDD'] = allMeanGrpVolRDD
    # get all attributes from dict and broadcast
    regDict['lambdaMotionGroups'] = lambdaMotion
    regDict['microMotionGainGroups'] = microMotionGain
    regDict['sigmaPixGroups'] = sigmaPix
    regDict['sigmaTimeGroups'] = sigmaTime
    regDict['structureSizeGroup'] = structureSize
    parIdx = regDict['parIdx']
    expandedTarget = regDict['expandedTarget']
    groupVolShiftsBC = sc.broadcast(regDict['grpVolShifts'])
    allMeanGrpVolRDD = regDict['allMeanGrpVolRDD']
    velPerPixXY = regDict['velPerPixXY']
    flyLines = regDict['flyLines']
    auPerPhoton = regDict['auPerPhoton']
    pixSizeXY = regDict['pixSizeXY']
    fieldTimes = regDict['fieldTimes']
    timeWeight = dict()
    for k in range(0, expandedTarget.shape[2]):
        timeWeight[k] = norm.pdf(np.absolute(fieldTimes[k] - fieldTimes) / (sigmaTime / microMotionGain)).reshape(-1,
                                                                                                                  1)
    timeWeightBC = sc.broadcast(timeWeight)
    bestTarget = copy.deepcopy(expandedTarget).astype(np.float64)
    for i in range(bestTarget.shape[2]):
        bestTarget[:, :, i] = inpaint_biharmonic(bestTarget[:, :, i], np.isnan(bestTarget[:, :, i]))
        # bestTarget[:, :, i] = white_tophat(bestTarget[:, :, i], structureSize)
    bestTargetBC = sc.broadcast(bestTarget)

    def getGroupWeightedPlanarDisplacement(kv):
        from numpy import unravel_index, argmax
        from skimage.feature import match_template, peak_local_max
        from skimage import measure
        from scipy.stats import norm
        from scipy.ndimage.interpolation import shift
        key, array1 = kv
        array1 = array1.astype(np.float64)
        bestTarget = bestTargetBC.value.astype(np.float64)

        initShift = -copy.deepcopy(groupVolShiftsBC.value[key[1], key[0], :-1])
        si2 = np.asarray(bestTarget.shape, dtype='float32')
        bestShift = np.zeros((int(si2[2]), 2))
        conWeight = np.zeros((int(si2[2]), 1))

        # initial plane by plane displacement
        for z in range(int(si2[2])):
            planeShift = copy.deepcopy(initShift)
            planeShift[0] = planeShift[0] - np.floor(flyLines[z] / 2).astype(int)
            grid0, grid1 = np.meshgrid(np.asarray(range(int(si2[1]))), np.asarray(range(int(si2[0]))))
            grid0 = (grid0 - np.floor(si2[1] / 2) + planeShift[1]).flatten()
            grid1 = (grid1 - np.floor(si2[0] / 2) + planeShift[0]).flatten()
            # one pixel shift is considered noise and ignored
            grid0sign = np.sign(grid0)
            grid1sign = np.sign(grid1)
            grid0 = (np.absolute(grid0) - 1.)
            grid1 = (np.absolute(grid1) - 1.)
            grid0[np.nonzero(grid0 < 0)[0]] = 0
            grid1[np.nonzero(grid1 < 0)[0]] = 0
            grid0 = grid0 * grid0sign
            grid1 = grid1 * grid1sign
            grid0 = grid0 * velPerPixXY
            grid1 = grid1 * velPerPixXY
            distGrid = (grid0 ** 2 + grid1 ** 2) ** 0.5
            distGrid = distGrid.reshape(np.round(si2[0:2]).astype(int))
            # convert to log-linear probability
            distGrid = np.exp(distGrid * (-np.asarray(1.0 / (lambdaMotion * microMotionGain), dtype='float64')))
            arrayPlane = copy.deepcopy(array1[:, :, z])
            arrayPlane = arrayPlane[flyLines[z]:, :]
            # arrayPlane = white_tophat(arrayPlane[flyLines[z]:, :], structureSize)

            targetPlane = bestTarget[:, :, z]

            # calculate shift probabilities
            nPhotonM = np.sum(arrayPlane.flatten()) / auPerPhoton
            c = match_template(targetPlane, arrayPlane, pad_input=True)

            # find max correlation
            s2 = c.shape
            maxC = np.nanmax(c.flatten())

            # for each shift in xcorr, convert to probability it is greater than actual peak Fisher correction
            maxC2 = 0.5 * np.log((1.0 + maxC) / (1.0 - maxC))
            c2 = 0.5 * np.log((1.0 + c) / (1.0 - c))

            # z-score, p-value
            zScore = (maxC2 - c2) / ((1.0 / (nPhotonM - 3) + 1.0 / (nPhotonM - 3)) ** 0.5)
            p = (1 - norm.cdf(zScore)) * 2

            # prevent weighted micro shifts by retaining only local maxima > 1um apart
            p2 = p * peak_local_max(p, min_distance=2000.0 / pixSizeXY, indices=False, threshold_rel=0,
                                    exclude_border=False).astype('float64')

            # weight against peaks that are improbably fast movements (distGrid)
            bestShiftPlane = unravel_index(argmax(p2 * distGrid), s2)
            bestShift[z, :] = -(bestShiftPlane - np.floor(np.asarray(p.shape) / 2))

            # Estimate the "gyration distance" around the selected peak using second moments of image
            # Moment Functions in Image Analysis: Theory and Applications
            # Mukundan and Ramakrishnan
            cm = measure.moments_central(p, bestShiftPlane[0], bestShiftPlane[1], order=2)
            xGyr = (cm[0, 2] / cm[0, 0]) ** 0.5
            yGyr = (cm[2, 0] / cm[0, 0]) ** 0.5
            eqGyr = (xGyr * yGyr / np.pi) ** 0.5
            conWeight[z] = (1 - norm.cdf(eqGyr / sigmaPix)) * 2

        # neighborhood displacement depending on registration certainty
        bestShift2 = copy.deepcopy(bestShift)
        for k2 in range(int(si2[2])):
            netWeight = timeWeightBC.value[k2] * conWeight
            bestShift2[k2, :] = bestShift[argmax(netWeight), :]
            bestShift2[k2, 0] = bestShift2[k2, 0] + np.floor(flyLines[k2] / 2).astype(int)

        # shift and NaN fly lines
        for z in range(int(array1.shape[2])):
            # if optimization:
            #     array1[:, :, z] = white_tophat(array1[:, :, z], structureSize)
            array1[array1 < 0] = 0
            array1[0:flyLines[z], :, z] = 0
            if useVolumeShifts:
                s = -initShift
                array1[:, :, z] = shift(array1[:, :, z], s, cval=float('nan'))
                bestShift2[z, :] = initShift
            else:
                s = -bestShift2[z, :]
                array1[:, :, z] = shift(array1[:, :, z], s, cval=float('nan'))
            nFlyLines = int(flyLines[z] + s[0])
            if nFlyLines > 0:
                array1[0:nFlyLines, :, z] = np.NAN

        # array1[array1 < 0] = np.NAN
        return key, (bestShift2, array1)

    # registration of groups to initial target
    regGrpVol = allMeanGrpVolRDD.map(getGroupWeightedPlanarDisplacement).collectAsMap()
    # pulling out the registered plains
    nPlane = regGrpVol[parIdx[0]][1].shape[2]
    planeGrpAlign = [None] * nPlane
    shifts = [None] * nPlane
    for j in range(int(nPlane)):
        nextIdx = np.asarray([i[0] == j for i in parIdx]).nonzero()[0]
        planeGrpAlign[j] = np.asarray([regGrpVol[parIdx[i]][1][:, :, j] for i in nextIdx])
        shifts[j] = np.asarray([regGrpVol[parIdx[i]][0][j, :] for i in nextIdx])
        planeGrpAlign[j] = planeGrpAlign[j].reshape(planeGrpAlign[j].shape[0], -1)

    regDict['planeGrpAlign'] = planeGrpAlign
    regDict['shiftsGroups'] = shifts

    if not optimization:
        # reorder shifts per time point
        groups = regDict['groups'].astype(int)
        shifts2 = np.stack(shifts).transpose(1, 0, 2)
        shiftsGroupsTime = np.zeros((groups.shape[0], groups.shape[1], 2), dtype=int)
        for timePoint, allPlaneGroups in enumerate(groups):
            for plane, group in enumerate(allPlaneGroups):
                shiftsGroupsTime[timePoint, plane, :] = shifts2[group, plane, :]

        regDict['shiftsGroupsTime'] = shiftsGroupsTime
        data = list()
        data.append(dict(y=np.max(regDict['shiftsGroupsTime'], axis=0)[:, 0], type='scatter', name='max x'))
        data.append(dict(y=np.max(regDict['shiftsGroupsTime'], axis=0)[:, 1], type='scatter', name='max y'))
        data.append(dict(y=np.mean(regDict['shiftsGroupsTime'], axis=0)[:, 0], type='scatter', name='mean x'))
        data.append(dict(y=np.mean(regDict['shiftsGroupsTime'], axis=0)[:, 1], type='scatter', name='mean y'))
        data.append(dict(y=np.std(regDict['shiftsGroupsTime'], axis=0)[:, 0], type='scatter', name='std x'))
        data.append(dict(y=np.std(regDict['shiftsGroupsTime'], axis=0)[:, 1], type='scatter', name='std y'))
        py.iplot(dict(data=data))
        m, s = divmod(time.time() - start, 60)
        logger.info('Time:  %02d:%02d' % (m, s))


@fallback
def getFinalGroupsSVD(sc, session, regDict, flyLinesFlag=False, finalClusterNum=45, finalIterKMeans=600,
                      finalNInit=25, finalTol=1e-8, finalKSVD=80, finalIterSVD=30, largeGroupMean=False, crop=None):
    """
    
    :param sc: Spark Context
    :param session: SpineSession object
    :param regDict: registration dict
    :param flyLinesFlag: if to exclude fly lines before shifting data in XY
    :param finalClusterNum: Number of clusters to form
    :param finalIterKMeans: Number of iteration in the KMeans
    :param finalNInit: Number if initialization in the KMeans int
    :param finalTol:Tolerance the KMeans
    :param finalKSVD: Number of SVD components to calculate
    :param finalIterSVD: Number of iteration in the SVD 
    :param largeGroupMean:
    :param crop: crop planes / lines with crossing cells. If None will not crop.
    :return: adds to regDict:
        1) finalGroups: cluster assignments
        2) finalGrpImg: the groups themselves
    """
    # get data from regDict
    start = time.time()
    regDict['finalClusterNum'] = finalClusterNum
    regDict['finalIterKMeans'] = finalIterKMeans
    regDict['finalNInit'] = finalNInit
    regDict['finalTol'] = finalTol
    regDict['finalKSVD'] = finalKSVD
    regDict['finalIterSVD'] = finalIterSVD
    regDict['finalCrop'] = crop
    indexNoAP = np.where(regDict['noAP'])[0]
    flyLines = regDict['flyLines']

    data = regDict['data']

    # broadcast shifts
    shiftsBC = sc.broadcast(regDict['shiftsGroupsTime'])

    def applyShifts_inner(kv):
        from scipy.ndimage.interpolation import shift as sp_shift
        from scipy.stats.mstats import zscore
        k, v = kv
        k = np.array(k[0]).astype(int)
        v = v.astype(np.float32)
        if extendedFlag:
            sz = v.shape
            out = np.zeros((sz[0] * 2, sz[1] * 2, sz[2])).astype('float32')
            out[:] = cVal
            if flyLinesFlag:
                for p in range(v.shape[2]):
                    v[:flyLines[p], :, p] = cVal
            out[sz[0] // 2:(sz[0] // 2 + sz[0]), sz[1] // 2:(sz[1] // 2 + sz[1]), :] = v
            v = out
        for p in range(v.shape[2]):
            shift = shiftsBC.value[k, p]
            plane = v[:, :, p]
            if flyLinesFlag and not extendedFlag:
                plane[0:flyLines[p], :] = -0.001
            plane = sp_shift(plane, -shift, cval=cVal, prefilter=False)
            if flyLinesFlag and not extendedFlag:
                plane[plane < 0] = cVal
            v[:, :, p] = plane
        if not np.isnan(cVal):
            v = zscore(v, axis=None)
        return v

    extendedFlag = False
    cVal = 0
    shiftedData = data.map(applyShifts_inner, with_keys=True)
    if crop is not None:
        shiftedData = shiftedData.map(lambda x: x[crop['ySlice'], crop['xSlice'], crop['zIndex']])
    shiftedData.cache()
    shiftedData.count()

    # train without time point with APs
    m, s = divmod(time.time() - start, 60)
    logger.info('Pre:  %02d:%02d' % (m, s))
    sys.stdout.flush()
    U, _, _ = getSVD(shiftedData, finalKSVD, getComponents=False, getS=False, normalization='mean')

    m, s = divmod(time.time() - start, 60)
    logger.info('SVD:  %02d:%02d' % (m, s))
    sys.stdout.flush()
    model = KMeans(n_clusters=finalClusterNum, init='k-means++', n_init=finalNInit, max_iter=finalIterKMeans,
                   tol=finalTol, precompute_distances='auto', verbose=0, random_state=None, copy_x=False, n_jobs=-1,
                   algorithm='auto').fit(U[indexNoAP, ...])

    clusters = model.predict(U)

    m, s = divmod(time.time() - start, 60)
    logger.info('Cluster:  %02d:%02d' % (m, s))

    # remove empty groups
    uIndex, loc = np.unique(clusters, return_inverse=True)
    clusters2 = np.zeros(clusters.shape)
    for i in range(len(uIndex)):
        clusters2[loc == i] = i
    clusters = clusters2

    # return data in expanded space
    extendedFlag = True
    cVal = float('NaN')
    returnData = data.map(applyShifts_inner, with_keys=True)
    sys.stdout.flush()

    try:
        # will fail due to Spark handling of data bigger then 2^32 see https://issues.apache.org/jira/browse/SPARK-6235
        finalGroup = nanMeanByIndex(returnData, clusters, largeGroupMean=largeGroupMean)
    except Exception:
        # doing it one cluster at a time
        logger.info('failed with nanMeanByIndex')
        returnData.cache()
        returnData.count()
        sz2 = returnData.shape
        finalGroup = np.zeros((len(uIndex), sz2[1], sz2[2], sz2[3]))
        for i, index in log_progress(enumerate(uIndex), every=1, size=len(uIndex), name='Groups'):
            pos = np.sort(np.where(clusters == index)[0])
            if len(pos) == 1:
                finalGroup[i, :, :, :] = returnData[pos, :, :, :].toarray()
            else:
                # logger.info('%d: %d, ' % (i, len(pos)))
                current = returnData[pos, :, :, :]
                current.cache()
                current.count()
                mean = np.nan_to_num(np.nanmean(current.toarray(), axis=0))
                MeanSz = mean.shape
                MeanVec = mean.reshape(1, MeanSz[0] * MeanSz[1] * MeanSz[2])
                CC = current.map(
                    lambda vol: np.corrcoef(np.nan_to_num(vol.reshape(1, MeanSz[0] * MeanSz[1] * MeanSz[2])),
                                            MeanVec)[0, 1]).toarray()
                CC = np.array(CC)
                length = current.shape[0]
                if length < 30:
                    points = np.array(length)
                elif length < 100:
                    points = np.round(length / 2)
                elif length < 200:
                    points = np.round(length / 3)
                elif length < 300:
                    points = np.round(length / 4)
                else:
                    points = np.round(length / 5)
                points = points.astype(int)
                ind = np.argpartition(CC, -points)[-points:]
                finalGroup[i, :, :, :] = np.nanmean(current[ind, :, :, :].toarray(), axis=0)
                current.uncache()
        returnData.uncache()

    regDict['finalGroups'] = clusters + 1
    regDict['finalGrpImg'] = finalGroup
    returnData.uncache()
    m, s = divmod(time.time() - start, 60)
    logger.info('Done:  %02d:%02d' % (m, s))
    sys.stdout.flush()

    # plot results
    ax1 = plt.subplot(1, 1, 1)
    index2 = np.unique(clusters)
    _, bins, _ = ax1.hist(clusters + 1, bins=len(index2))
    ax1.set_ylabel('Timepoints', color='b')
    ax1.set_xlabel('Cluster#')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    binCenters = 0.5 * (bins[1:] + bins[:-1])
    ax2 = ax1.twinx()
    Mean = np.nanmean(finalGroup, axis=(1, 2, 3))
    Std = np.nanstd(finalGroup, axis=(1, 2, 3))
    ax2.errorbar(binCenters, Mean, Std, ecolor='r', fmt='ro')
    plt.grid(False)
    ax2.set_ylabel(u'Intensity (±SD)', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    plt.show()
    return sc


def estZShiftsGrouping(sc, regDict, nZShifts=5, zFWHM=2.5, nanThreshold=200, crop=None, preFlag=True,
                       redoBestCrop=False, do_plot=True):
    """ estimate Z shifts and make nZshifts new groups from postMeanGrpVol

    :param sc: Spark Context
    :param regDict: needs: postMeanGrpVol, Index, nCutLines, pixSizeXY, finalGroups, newShifts
    :param nZShifts: final number of Z shifts to coalesce data into
    :param zFWHM: FWHM of the z-axis PSF
    :param nanThreshold
    :param crop: crop dict from getCrop, if None will crop 4 lines in X and Y
    :param preFlag: flag for doing tophat and z score on the groups
    :param redoBestCrop: redo the best crop step, if false will take it from regDict
    :return: adds: groupZ, grpZIdx, grpZPos, finalShifts
    """
    start = time.time()

    regDict['nZShifts'] = nZShifts
    regDict['nanThresholdZShiftsGrouping'] = nanThreshold
    regDict['zFWHM'] = zFWHM
    groups = copy.deepcopy(regDict['finalGrpImg'])
    finalGroups = copy.deepcopy(regDict['finalGroups'])

    ngrp = np.array([np.sum(finalGroups == grp) for grp in np.unique(finalGroups)]).astype('float32')
    weightedGroups = np.zeros(groups.shape)
    for i in range(groups.shape[0]):
        weightedGroups[i, :, :, :] = groups[i, :, :, :] * ngrp[i]
    # get the best crop based on number of NaNs
    origSz = groups.shape
    logger.info('Original size: %s' % (origSz,))
    groups = groups.transpose(1, 2, 3, 0).reshape(origSz[1], origSz[2], origSz[0] * origSz[3])
    maxShifts = regDict['shiftsGroupsTime'].max(axis=(0, 1))
    if maxShifts[0] > origSz[1] / 4:
        maxShifts[0] = origSz[1] / 4
    if maxShifts[1] > origSz[2] / 4:
        maxShifts[1] = origSz[2] / 4
    if maxShifts[0] < 4:
        maxShifts[0] = 4
    if maxShifts[1] < 4:
        maxShifts[1] = 4

    # check if crop is already made and used cached version
    if 'groupsInit' not in regDict.keys() or redoBestCrop:
        logger.info('Generating best crop')
        returned, groups = getBestCrop(sc, groups, maxShifts=maxShifts, nanThreshold=nanThreshold * origSz[0],
                                       inPaint=True, checkSize=False)
        regDict['returnedInit'] = copy.deepcopy(returned)
        regDict['groupsInit'] = copy.deepcopy(groups)
    else:
        logger.info('Using cached best crop')
        returned = copy.deepcopy(regDict['returnedInit'])
        groups = copy.deepcopy(regDict['groupsInit'])

    expSz = groups.shape
    groups = groups.reshape(expSz[0], expSz[1], origSz[3], origSz[0]).transpose(3, 0, 1, 2)

    # crop cell body and some lines
    if crop is None:
        lines = 4
        logger.info('No crop given, using %d lines' % lines)
        groups = groups[:, lines:-lines, lines:-lines, :]
        regDict['estZCropLines'] = lines
        regDict['estZCropDict'] = None
    else:
        groups = groups[:, crop['ySlice'], crop['xSlice'], crop['zIndex']]
        regDict['estZCropLines'] = None
        regDict['estZCropDict'] = crop
    m, s = divmod(time.time() - start, 60)
    logger.info('group shape after crop: %s, time: %02d:%02d' % (groups.shape, m, s))

    # calculate the correlation dissimilarity between all groups
    from scipy.stats.mstats import zscore
    groups2 = copy.deepcopy(groups)
    if preFlag:
        for i in range(groups.shape[0]):
            for j in range(groups.shape[3]):
                groups2[i, :, :, j] = white_tophat(groups2[i, :, :, j], 6)
            groups2[i, :, :, :] = zscore(groups2[i, :, :, :], axis=None)

    flatGroups = groups2.reshape(groups.shape[0], -1)
    R = np.corrcoef(np.nan_to_num(flatGroups))
    R = 1 - R

    # order the groups based on similarity using a traveling salesman problem solver
    bestPath = np.asarray(solve_tsp(R))
    regDict['bestPath'] = copy.deepcopy(bestPath)
    # Coalesce groups into nZShifts final zShift groups to minimize activity-based dissimilarity

    grpsPerShift = np.floor(bestPath.shape[0] / nZShifts).astype(int)
    si = weightedGroups.shape

    logger.info('weightedGroup shape: %s' % (si,))
    groupZ = np.zeros((nZShifts, si[1], si[2], si[3]))
    grpZIdx = [None] * nZShifts
    for i in range(nZShifts - 1):
        grpZIdx[i] = bestPath[(i * grpsPerShift):((i * grpsPerShift) + grpsPerShift)]
        groupZ[i, :, :, :] = np.nanmean(weightedGroups[grpZIdx[i], :, :, :], axis=0) * len(grpZIdx[i]) / np.sum(
            ngrp[grpZIdx[i]])
    grpZIdx[nZShifts - 1] = bestPath[((nZShifts - 1) * grpsPerShift):]
    groupZ[nZShifts - 1, :, :, :] = np.nanmean(weightedGroups[bestPath[((i + 1) * grpsPerShift):], :, :, :],
                                               axis=0) * len(grpZIdx[i]) / np.sum(
        ngrp[bestPath[((i + 1) * grpsPerShift):]])

    regDict['groupZfull'] = copy.deepcopy(groupZ)
    # crop groupZ as groups
    shiftIndex = np.argmax(np.array([s[1] for s in returned]))
    borderFinal = returned[shiftIndex][0]
    groupZ = groupZ[:, borderFinal[0]:-borderFinal[1], borderFinal[2]:-borderFinal[3], :]

    groupZCC = copy.deepcopy(groupZ)
    if crop is None:
        lines = 12
        groupZCC = groupZCC[:, lines:-lines, lines:-lines, :]
    else:
        groupZCC = groupZCC[:, crop['ySlice'], crop['xSlice'], crop['zIndex']]
    origVols = groupZCC
    # based on XY-shift correlations estimate Zshift between furthest groups
    for i in range(groupZCC.shape[0]):
        for k in range(groupZCC.shape[3]):
            origVols[i, :, :, k] = inpaint_biharmonic(groupZCC[i, :, :, k], np.isnan(groupZCC[i, :, :, k]))

    # calculate dissimilarity
    flatOrigVols = origVols.reshape(origVols.shape[0], -1)
    dissim = 1 - np.corrcoef(np.nan_to_num(flatOrigVols))

    # filter XY to the larger Z-PSF
    pixelSizeXY = regDict['pixSizeXY'] / float(1000)
    filterFWHM = zFWHM / pixelSizeXY
    filterSigma = filterFWHM / 2.355
    filtVols = np.zeros(origVols.shape)
    for i in range(0, origVols.shape[0]):
        for j in range(0, origVols.shape[3]):
            filtVols[i, :, :, j] = gaussian_filter(origVols[i, :, :, j], filterSigma, mode='wrap')
    # calculate autocorrelation
    ccXY = np.zeros(filtVols.shape[0:3])
    midPlane = np.floor(filtVols.shape[3] / 2).astype(int)
    for i in range(0, filtVols.shape[0]):
        ccXY[i, :, :] = 1 - match_template(filtVols[i, :, :, :], filtVols[i, :, :, :],
                                           pad_input=True, mode='wrap')[:, :, midPlane]
    # estimate shift to corr curve
    si = ccXY[0, :, :].shape
    grid0, grid1 = np.meshgrid(np.asarray(range(0, si[1])), np.asarray(range(0, si[0])))
    grid0 = (grid0 - np.floor(si[1] / 2)).flatten() * pixelSizeXY
    grid1 = (grid1 - np.floor(si[0] / 2)).flatten() * pixelSizeXY
    distGrid = (grid0 ** 2 + grid1 ** 2) ** 0.5
    distGrid = (distGrid.reshape(-1, 1) * np.ones((1, ccXY.shape[0]))).T.flatten()
    ccXY = ccXY.flatten()
    uniqueDist = np.unique(distGrid.flatten())
    ccDist = np.asarray([np.median(ccXY[distGrid == dist]) for dist in uniqueDist])
    idx = np.argsort(ccDist)
    f = PchipInterpolator(ccDist[idx], uniqueDist[idx])
    maxZShift = f(dissim[0, -1]) * 2

    # assign Z positions to each group assuming even spacing across groupZ
    zSpacing = maxZShift / (nZShifts - 1)
    grpZPos = (np.asarray(range(0, nZShifts)) - ((nZShifts - 1) / float(2))) * zSpacing

    # newShifts = regDict['grpVolShifts'][:,:,:-1]
    zshifts = np.ones(len(finalGroups))
    zshiftsPos = np.ones(len(finalGroups))
    for i in range(0, nZShifts):
        for j in grpZIdx[i]:
            zshiftsPos[finalGroups == (j + 1)] = grpZPos[i]
            zshifts[finalGroups == (j + 1)] = i

    # inpaint expandedZTargets
    origSz = groupZ.shape
    logger.info('origSz: %s' % (origSz,))
    groupZ = groupZ.transpose(1, 2, 3, 0).reshape(origSz[1], origSz[2], origSz[0] * origSz[3])

    def inPaintPar(volume):
        mask = np.isnan(volume)
        if np.sum(mask):
            try:
                return inpaint_biharmonic(volume, mask)
            except ValueError:
                return np.nan_to_num(volume)
        else:
            return volume

    groupZ = td.images.fromarray(groupZ, engine=sc).map(inPaintPar).toarray()
    expSz = groupZ.shape
    logger.info('expSz shape: %s' % (expSz,))
    groupZ = groupZ.reshape(expSz[0], expSz[1], origSz[3], origSz[0]).transpose(3, 0, 1, 2)

    regDict['grpZIdxGroup'] = grpZIdx
    regDict['grpZPosGroup'] = grpZPos
    regDict['zshiftsGroup'] = zshifts.astype(int)
    regDict['zshiftsPosGroup'] = zshiftsPos
    regDict['expandedZTargets'] = groupZ
    regDict['returnedEstZ'] = returned
    if do_plot:
        z, c = np.unique(regDict['zshiftsGroup'], return_counts=True)
        fig = tls.make_subplots(rows=1, cols=3, specs=[[{'colspan': 2}, None, {}]], print_grid=False)
        fig.append_trace(go.Scatter({'y': regDict['zshiftsGroup']}, ), 1, 1)
        x = list(map(str, regDict['grpZPosGroup']))
        fig.append_trace(go.Bar({'x': x, 'y': c, }, ), 1, 3)
        fig['layout']['xaxis1'].update(title='Time points')
        fig['layout']['yaxis1'].update(title='Z position (um)')
        fig['layout']['xaxis2'].update(title='Z Position')
        fig['layout']['yaxis2'].update(title='# Time points')
        fig['layout'].update(showlegend=False)
        py.iplot(fig)


def getShiftOffset(regDict):
    """

    :param regDict:
    :return:
    """
    shiftIndex = np.argmax(np.array([s[1] for s in regDict['returnedEstZ']]))
    borderFinal = np.array(regDict['returnedEstZ'][shiftIndex][0]).astype(float)
    dims = regDict['dims'][:2]
    mid = (np.array(dims) / 2).astype(float)
    left = (mid[0] - borderFinal[0])
    right = mid[0] - borderFinal[1]
    up = mid[1] - borderFinal[2]
    down = mid[1] - borderFinal[3]
    offset = np.array([int((left - right) / 2), int((up - down) / 2)])
    regDict['shiftsGroupsTime2'] = copy.deepcopy(regDict['shiftsGroupsTime']) + offset


def runVolReg(sc, regDict, TForm1, numTimepoints=10000, percentile=33):
    """

    :param sc: Spark Context
    :param TForm1: Transformation dict for embedding to 3d space
    :param regDict: registration dict
    :param numTimepoints: number of time points to use in optimization
    :param percentile: cutoff for the distribution os the optimization metric
    :return:
    """
    WeightedShiftsOptimization(sc, regDict, numTimepoints=numTimepoints, usePerPlanePrior=True, percentile=percentile)
    shiftList = ['newShifts', 'grpVolShiftsTime', 'shiftsGroupsTime']
    for shift in shiftList:
        logger.info(shift)
        shifts = regDict[shift]
        regTemp = applyShifts(sc, regDict, optimization=False, shifts=shifts)
        meanTemp = regTemp.map(lambda x: np.nan_to_num(x)).mean().toarray()
        embedTemp = getExampleVol(sc, meanTemp, TFormDict=TForm1, project=True)
        writeTiff(regDict['fullPath'], embedTemp, shift + '_regData')
    getAlignedGroups(regDict)


def WeightedShiftsOptimization(sc, regDict, numTimepoints=10000, usePerPlanePrior=False, percentile=33):
    """ optimization of time points alignment in x y by minimizing fano factor and penalizing non-linearly the number
        of NaNs, uses numTimepoints random subset of time points
    :param sc: Spark context
    :param regDict: registration dict
    :param numTimepoints: the number of time points to use
    :param usePerPlanePrior: if True will use the priors (False)
    :param percentile: threshold for optimization
    :return:
    """
    start = time.time()
    getShiftOffset(regDict)
    data = regDict['data']
    sz = data.shape[0]
    index = np.array(sample(range(sz), numTimepoints))
    regDict['dataMidIndex'] = index
    dataMid = data[index, :, :, :]
    dataMid = balanced_repartition(dataMid.tordd(), sc.defaultParallelism * 3)
    dataMid = td.images.fromrdd(dataMid.sortByKey())
    dataMid.cache()
    regDict['dataMid'] = dataMid
    m, s = divmod(time.time() - start, 60)
    logger.info('Got data for optimization with %d time points, %02d:%02d' % (numTimepoints, m, s))
    # get volume shift values
    array = applyShifts(sc, regDict, shifts=regDict['shiftsGroupsTime'][index, ...], optimization=True).toarray()
    Mean = np.nanmean(array, axis=0)
    t = np.nanpercentile(Mean, 40)
    mask = array < t
    sz2 = array.shape
    regDict['Nans_volume'] = np.sum(np.isnan(array)).astype(float) / np.prod(sz2) * 100
    array[mask] = np.nan
    mean = np.nanmean(array, axis=0)
    var = np.nanvar(array, axis=0)
    regDict['fano_volume'] = np.nanpercentile(var / mean, percentile)
    m, s = divmod(time.time() - start, 60)
    logger.info('Volume shifts fano: %.2f, Nans: %.2f, %02d:%02d' %
                (regDict['fano_volume'], regDict['Nans_volume'], m, s))

    if 'limitedShiftsLPBC' not in regDict.keys():
        groupedTargetBC = sc.broadcast(regDict['expandedZTargets'])
        regDict['groupedTargetBC'] = groupedTargetBC
        zshiftsBC = sc.broadcast(regDict['zshiftsGroup'])
        regDict['zshiftsBC'] = zshiftsBC
        regDict['shiftsGroupsTimeBC'] = sc.broadcast(regDict['shiftsGroupsTime2'])

    optimizationResult = differential_evolution(optimizeFano, bounds=[(0.5, 100), (3, 40), (10, 1000)],
                                                args=(sc, regDict, percentile, usePerPlanePrior, start),
                                                maxiter=5, popsize=3, mutation=(0.5, 1),
                                                disp=True,
                                                polish=False)
    regDict['optimizationResultTimepoints'] = optimizationResult
    regDict['optimizationResultTimepointsPercentile'] = percentile
    lambdaMotion = optimizationResult.x[0] / 100.0
    sigmaTime = optimizationResult.x[1]
    sigmaPix = optimizationResult.x[2] / 100.0
    getWeightedShifts(sc=sc, regDict=regDict, usePerPlanePrior=usePerPlanePrior, lambdaMotion=lambdaMotion,
                      sigmaPix=sigmaPix, sigmaTime=sigmaTime)
    regData = applyShifts(sc=sc, regDict=regDict)
    regDict['regData'] = regData
    regDict['regData'].cache()
    regDict['regData'].count()


def optimizeFano(params, sc, regDict, percentile, usePerPlanePrior, start):
    """ test params with getWeightedShifts by calculating the % change in # of NaNs and fano factor

    :param params: tuple of (lambda_motion, sigma_time, sigma_pix)
    :param sc: Spark Context
    :param regDict: registration dict
    :param percentile: to use for fano factor
    :param usePerPlanePrior: use priors in the registration function
    :param start: start time of the optimization function
    :return:
    """
    lambda_motion = float(params[0]) / 100.0
    sigma_time = float(params[1])
    sigma_pix = float(params[2]) / 100.0
    getWeightedShifts(sc, regDict, optimization=True, usePerPlanePrior=usePerPlanePrior,
                      lambdaMotion=lambda_motion, sigmaPix=sigma_pix, sigmaTime=sigma_time)
    after_shifts = applyShifts(sc, regDict, optimization=True)
    after_shifts.cache()
    after_shifts.count()
    Mean = after_shifts.map_as_series(lambda x: np.array(np.nanmean(x), ndmin=1)).toarray()
    sz2 = after_shifts.shape
    Nans = np.sum(after_shifts.map(lambda x: np.sum(np.isnan(x))).sum().toarray()).astype(float) / np.prod(sz2) * 100
    t = np.nanpercentile(Mean, 40)

    def apply_mask(x):
        mask = x < t
        x[mask] = np.nan
        return x

    mean = after_shifts.map(apply_mask).map_as_series(lambda x: np.array(np.nanmean(x), ndmin=1)).toarray()
    var = after_shifts.map(apply_mask).map_as_series(lambda x: np.array(np.nanvar(x), ndmin=1)).toarray()
    fano = np.nanpercentile(var / mean, percentile)

    assert fano > 0.0
    fano_rel = fano / regDict['fano_volume'] * 100 - 100
    Nans_rel = Nans - regDict['Nans_volume']
    result = fano_rel + Nans_rel
    after_shifts.uncache()
    m, s = divmod(time.time() - start, 60)
    logger.info(
        'L:% 6.3f,T:% 4.1f,P:% 4.2f,N:%.2f,f:%.2f,r:%.3f, %02d:%02d' %
        (lambda_motion, sigma_time, sigma_pix, Nans_rel, fano_rel, result, m, s))
    return result


def getWeightedShifts(sc, regDict, optimization=False, usePerPlanePrior=False, lambdaMotion=2.0, microMotionGain=1.0,
                      sigmaPix=3.0, sigmaTime=36.0):
    """ X Y align time point while considering certainty, prior and neighboring planes

    :param sc: Spark Context
    :param regDict: registration dict
    :param optimization: if True called from optimization function WeightedShiftsOptimization
    :param usePerPlanePrior: if True will use priors instead of calculated shifts
    :param lambdaMotion: How big of a motion do we want to consider for current vs. prior
    :param microMotionGain: currently not used (set to 1.0)
    :param sigmaPix: How big of a motion do we want to consider for current vs. previous plane
    :param sigmaTime: how much do we weight the previous shifts in time
    :return: newShifts
    """
    if not optimization:
        # limitedShiftsLPBC = sc.broadcast(regDict['limitedShiftsLP'])
        # regDict['limitedShiftsLPBC'] = limitedShiftsLPBC
        groupedTargetBC = sc.broadcast(regDict['expandedZTargets'])
        regDict['groupedTargetBC'] = groupedTargetBC
        zshiftsBC = sc.broadcast(regDict['zshiftsGroup'])
        regDict['zshiftsBC'] = zshiftsBC
        shiftsGroupsTimeBC = sc.broadcast(regDict['shiftsGroupsTime2'])
        regDict['shiftsGroupsTimeBC'] = shiftsGroupsTimeBC
    else:
        # limitedShiftsLPBC = regDict['limitedShiftsLPBC']
        groupedTargetBC = regDict['groupedTargetBC']
        zshiftsBC = regDict['zshiftsBC']
        shiftsGroupsTimeBC = regDict['shiftsGroupsTimeBC']

    dims = regDict['expandedZTargets'][0, :, :, :].shape
    velPerPixXY = regDict['velPerPixXY']
    regDict['lambdaMotionTimepoints'] = lambdaMotion
    regDict['microMotionGainTimepoints'] = microMotionGain
    regDict['sigmaPixTimepoints'] = sigmaPix
    regDict['sigmaTimeTimepoints'] = sigmaTime
    medianFly = np.median(regDict['flyLines'])
    auPerPhoton = regDict['auPerPhoton']
    pixSizeXY = regDict['pixSizeXY']
    flyLines = copy.deepcopy(regDict['flyLines'])
    fieldMaskBool = regDict['fieldMaskBool']
    fieldDur = regDict['fieldDur']
    optimizationIndex = regDict['dataMidIndex']
    # pre calculate time and
    fieldTimes = np.asarray(range(len(fieldMaskBool))) * fieldDur
    fieldTimes = fieldTimes[fieldMaskBool]
    timeWeight = dict()
    for k in range(regDict['dims'][2]):
        timeWeight[k] = norm.pdf(np.absolute(fieldTimes[k] - fieldTimes) / (sigmaTime / microMotionGain)).reshape(-1, 1)

    # offset per plane
    offset = (medianFly / 2).astype(int)
    medianFly = int(medianFly)
    nanLines = flyLines - medianFly
    si = np.asarray(dims, dtype='float32')

    def getWeightedPlanarDisplacement(kv):
        from numpy import unravel_index, argmax
        from skimage.feature import match_template, peak_local_max
        from skimage import measure
        from scipy.stats import norm

        key, array1 = kv
        key = np.array(key[0]).astype(int)
        if optimization:
            key = (optimizationIndex[key],)
        array1 = array1.astype(dtype='float64')
        array1 = array1[medianFly:, :, :]
        for i in range(array1.shape[2]):
            if nanLines[i] > 0:
                array1[:nanLines[i], :, i] = 0.0
                # array1[:, :, i] = inpaint_biharmonic(array1[:, :, i], np.isnan(array1[:, :, i]))

        bestShift = np.zeros((si[2].astype(int), 2))
        conWeight = np.zeros((si[2].astype(int), 1))

        # get z Target
        bestTarget = groupedTargetBC.value[zshiftsBC.value[key[0]], :, :, :].squeeze()

        # initial plane by plane displacement
        for i in range(int(si[2])):
            if usePerPlanePrior:
                initShift = shiftsGroupsTimeBC.value[key, i]
            # else:
            #     initShift = limitedShiftsLPBC.value[key[0]]
            initShift[0] = initShift[0] - offset
            a = range(int(si[1]))
            b = range(int(si[0]))
            grid0, grid1 = np.meshgrid(np.asarray(a), np.asarray(b))
            grid0 = (grid0 - np.floor(si[1] / 2) + initShift[1]).flatten()
            grid1 = (grid1 - np.floor(si[0] / 2) + initShift[0]).flatten()
            # one pixel shift is considered noise and ignored
            grid0sign = np.sign(grid0)
            grid1sign = np.sign(grid1)
            grid0 = (np.absolute(grid0) - 1.)
            grid1 = (np.absolute(grid1) - 1.)
            grid0[np.nonzero(grid0 < 0)[0]] = 0
            grid1[np.nonzero(grid1 < 0)[0]] = 0
            grid0 = grid0 * grid0sign
            grid1 = grid1 * grid1sign
            grid0 = grid0 * velPerPixXY
            grid1 = grid1 * velPerPixXY
            distGrid = (grid0 ** 2 + grid1 ** 2) ** 0.5
            distGrid = distGrid.reshape(np.array(si[0:2]).astype(int))
            # convert to log-linear probability
            distGrid = np.exp(distGrid * (-np.array(1.0 / (lambdaMotion * microMotionGain), dtype='float64')))

            arrayPlane = array1[:, :, i]
            # arrayPlane[0:flyLines[i], :] = 0
            targetPlane = bestTarget[:, :, i]

            # calculate shift probabilities
            nPhotonM = np.sum(arrayPlane.flatten()) / auPerPhoton
            c = match_template(targetPlane, arrayPlane, pad_input=True)

            # find max correlation
            s = c.shape
            maxCorr = np.nanmax(c.flatten())

            # for each shift in xcorr, convert to probability it is greater than actual peak
            # Fisher correction
            maxCorr2 = 0.5 * np.log((1.0 + maxCorr) / (1.0 - maxCorr))
            c2 = 0.5 * np.log((1.0 + c) / (1.0 - c))

            # z-score
            zScore = (maxCorr2 - c2) / ((1.0 / (nPhotonM - 3) + 1.0 / (nPhotonM - 3)) ** 0.5)

            # p-value
            p = (1 - norm.cdf(zScore)) * 2

            # prevent weighted microshifts by retaining only local maxima > 2um apart
            p2 = p * peak_local_max(p, min_distance=float(2000) / pixSizeXY, indices=False, threshold_rel=0,
                                    exclude_border=False).astype('float64')

            # weight against peaks that are improbably fast movements (distGrid)
            bestShiftPlane = unravel_index(argmax(p2 * distGrid), s)
            bestShift[i, :] = -(bestShiftPlane - np.floor(np.asarray(p.shape) / 2))

            # Estimate the "gyration distance" around the selected peak using second moments of image
            # Moment Functions in Image Analysis: Theory and Applications
            # Mukundan and Ramakrishnan
            cm = measure.moments_central(p, bestShiftPlane[0], bestShiftPlane[1], order=2)
            xGyr = (cm[0, 2] / cm[0, 0]) ** 0.5
            yGyr = (cm[2, 0] / cm[0, 0]) ** 0.5
            eqGyr = (xGyr * yGyr / np.pi) ** 0.5
            conWeight[i] = (1 - norm.cdf(eqGyr / sigmaPix)) * 2

        # neighborhood displacement depending on registration certainty
        bestShift2 = bestShift
        for i in range(int(si[2])):
            netWeight = timeWeight[i] * conWeight
            bestShift2[i, :] = bestShift[argmax(netWeight), :]

        bestShift2[:, 0] = bestShift2[:, 0] + offset
        return bestShift2

    if optimization:
        newShifts = regDict['dataMid'].map(getWeightedPlanarDisplacement, with_keys=True, dtype=np.float64,
                                           value_shape=(regDict['dims'][2], 2)).toarray()
    else:
        newShifts = regDict['data'].map(getWeightedPlanarDisplacement, with_keys=True, dtype=np.float64,
                                        value_shape=(regDict['dims'][2], 2)).toarray()
    regDict['newShifts'] = newShifts
    if not optimization:
        plt.figure(figsize=(10, 5))
        meanX = np.mean(newShifts[:, :, 0], axis=1)
        meanY = np.mean(newShifts[:, :, 1], axis=1)
        stdX = np.std(newShifts[:, :, 0], axis=1)
        stdY = np.std(newShifts[:, :, 1], axis=1)
        plt.subplot(2, 2, 1)
        plt.plot(meanX, 'b')
        plt.ylabel('mean x')
        plt.subplot(2, 2, 2)
        plt.plot(meanY, 'b')
        plt.ylabel('mean y')
        plt.subplot(2, 2, 3)
        plt.plot(stdX, 'r')
        plt.ylabel('std x')
        plt.subplot(2, 2, 4)
        plt.plot(stdY, 'r')
        plt.ylabel('std y')


# regDict['limitedShiftsLP'][:,:2]
def applyShifts(sc, regDict, optimization=False, shifts=None, useGroupShifts=False):
    """ shift the images object according to the new shifts

    :param sc: SparkContext
    :param regDict: needs newShifts, flyLines
    :param optimization: if True was called from WeightedShiftsOptimization: will take the dataMid and will also top hat
                         the result
    :param shifts: which shifts to use
    :param useGroupShifts: use limitedShiftsLP - the prior
    :return: registered Images object
    """
    if shifts is None:
        shifts = regDict['newShifts']
    flyLines = regDict['flyLines']
    AllShiftsArrayBC = sc.broadcast(shifts)

    def XYVolumePlanarAlign(kv):
        key, im = kv
        key = np.array(key[0]).astype(int)
        from scipy.ndimage.interpolation import shift
        im = im.astype('float32')
        for z in range(0, im.shape[2]):
            # if optimization:
            #     im[:, :, z] = white_tophat(im[:, :, z], 3)
            #     im[im < 0] = 0
            im[0:flyLines[z], :, z] = 0
            if useGroupShifts:
                s = -AllShiftsArrayBC.value[key, :]
                im[:, :, z] = shift(im[:, :, z], s, cval=float('nan'))
            else:
                s = -AllShiftsArrayBC.value[key, z, :]
                im[:, :, z] = shift(im[:, :, z], s, cval=float('nan'))
                nFlyLines = int(flyLines[z] + s[0])
                if nFlyLines > 0:
                    im[0:nFlyLines, :, z] = np.NAN

        return im.astype(dtype='float32')

    if optimization:
        return regDict['dataMid'].map(XYVolumePlanarAlign, dtype=np.float32, value_shape=regDict['dims'],
                                      with_keys=True)
    else:
        return regDict['data'].map(XYVolumePlanarAlign, dtype=np.float32, value_shape=regDict['dims'], with_keys=True)


def getAlignedGroups(regDict):
    """ Get groups after x y alignment with (postMeanGrpVolNo) or without APs (postMeanGrpVolNoAP)

    :param regDict: registration dict
    :return: adds: postMeanGrpVol, postMeanGrpVolNoAP
    """
    start = time.time()
    finalGroups = copy.deepcopy(regDict['finalGroups'])

    # Same group averaging, but replace time-points during putative APs with last preceding non-AP time-point
    expPutAP = signal.convolve(regDict['putAP'], [1, 1, 1], mode='valid') > 0
    notAP = np.logical_not(expPutAP).nonzero()[0]
    reIdx = np.array(range(len(finalGroups)))
    for i in range(len(finalGroups)):
        newIdx = (notAP <= i).nonzero()[0]
        if np.any(newIdx):
            reIdx[i] = notAP[newIdx[-1]]

    postMeanGrpVol = nanMeanByIndex(regDict['regData'], finalGroups)
    postMeanGrpVolNoAP = nanMeanByIndex(regDict['regData'], finalGroups[reIdx])

    regDict['postMeanGrpVol'] = postMeanGrpVol
    regDict['postMeanGrpVolNoAP'] = postMeanGrpVolNoAP
    m, s = divmod(time.time() - start, 60)
    logger.info('getAlignedGroups:  %02d:%02d' % (m, s))


def estZShifts(regDict, nZShifts, zFWHM, crop, noAP=True, shifts='newShifts'):
    """ estimate Z shifts and make nZshifts new groups from postMeanGrpVol

    :param regDict: needs: postMeanGrpVol, Index, nCutLines, pixSizeXY, finalGroups, newShifts
    :param nZShifts: final number of Z shifts to coalesce data into
    :param zFWHM: FWHM of the z-axis PSF
    :param crop: a cropping to apply before estimating the z positions
    :param noAP: exclude time points with APs
    :param shifts: name of shifts to take: newShifts - after optimization, shiftsGroupsTime - before
    :return: adds: groupZ, grpZIdx, grpZPos, finalShifts
    """

    regDict['nZShifts'] = nZShifts
    regDict['zFWHM'] = zFWHM
    if noAP:
        groups = copy.deepcopy(regDict['postMeanGrpVolNoAP'])
    else:
        groups = copy.deepcopy(regDict['postMeanGrpVol'])
    finalGroups = copy.deepcopy(regDict['finalGroups'])

    ngrp = np.array([np.sum(finalGroups == grp) for grp in np.unique(finalGroups)])
    weightedGroups = np.zeros(groups.shape)
    for i in range(groups.shape[0]):
        weightedGroups[i, :, :, :] = groups[i, :, :, :] * ngrp[i]

    groups = groups[:, crop['ySlice'], crop['xSlice'], crop['zIndex']]

    logger.info('group shape:' + str(groups.shape))
    # calculate the correlation dissimilarity between all groups
    flatGroups = groups.reshape(groups.shape[0], -1)
    R = np.corrcoef(np.nan_to_num(flatGroups))
    R = 1 - R

    # order the groups based on similarity using a traveling salesman problem solver
    bestPath = np.asarray(solve_tsp(R))

    # Coalesce groups into nZShifts final zShift groups to minimize activity-based dissimilarity
    grpsPerShift = np.floor(bestPath.shape[0] / nZShifts).astype(int)
    si = weightedGroups.shape
    groupZ = np.zeros((nZShifts, si[1], si[2], si[3]))
    grpZIdx = [None] * nZShifts
    for i in range(0, nZShifts - 1):
        grpZIdx[i] = bestPath[(i * grpsPerShift):((i * grpsPerShift) + grpsPerShift)]
        groupZ[i, :, :, :] = np.nansum(weightedGroups[grpZIdx[i], :, :, :], axis=0) / np.sum(ngrp[grpZIdx[i]])
    grpZIdx[nZShifts - 1] = bestPath[((nZShifts - 1) * grpsPerShift):]
    groupZ[nZShifts - 1, :, :, :] = np.nansum(weightedGroups[bestPath[(i * grpsPerShift):], :, :, :],
                                              axis=0) / np.sum(ngrp[bestPath[(i * grpsPerShift):]])

    # add registration of groupZ from center to edges

    groupZCC = copy.deepcopy(groupZ)
    if crop is None:
        lines = 12
        groupZCC = groupZCC[:, lines:-lines, lines:-lines, :]
    else:
        groupZCC = groupZCC[:, crop['ySlice'], crop['xSlice'], crop['zIndex']]
    # based on XY-shift correlations estimate Zshift between furthest groups
    origVols = np.nan_to_num(groupZCC)

    # calculate dissimilarity
    flatOrigVols = origVols.reshape(origVols.shape[0], -1)
    dissim = 1 - np.corrcoef(np.nan_to_num(flatOrigVols))

    # filter XY to the larger Z-PSF
    pixelSizeXY = regDict['pixSizeXY'] / float(1000)
    filterFWHM = zFWHM / pixelSizeXY
    filterSigma = filterFWHM / 2.355
    filtVols = np.zeros(origVols.shape)
    for i in range(0, origVols.shape[0]):
        for j in range(0, origVols.shape[3]):
            filtVols[i, :, :, j] = gaussian_filter(origVols[i, :, :, j], filterSigma, mode='wrap')
    # calculate autocorrelation
    ccXY = np.zeros(filtVols.shape[0:3])
    midPlane = np.floor(filtVols.shape[3] / 2).astype(int)
    for i in range(0, filtVols.shape[0]):
        ccXY[i, :, :] = 1 - match_template(filtVols[i, :, :, :], filtVols[i, :, :, :],
                                           pad_input=True, mode='wrap')[:, :, midPlane]
    # estimate shift to corr curve
    si = ccXY[0, :, :].shape
    grid0, grid1 = np.meshgrid(np.asarray(range(0, si[1])), np.asarray(range(0, si[0])))
    grid0 = (grid0 - np.floor(si[1] / 2)).flatten() * pixelSizeXY
    grid1 = (grid1 - np.floor(si[0] / 2)).flatten() * pixelSizeXY
    distGrid = (grid0 ** 2 + grid1 ** 2) ** 0.5
    distGrid = (distGrid.reshape(-1, 1) * np.ones((1, ccXY.shape[0]))).T.flatten()
    ccXY = ccXY.flatten()
    uniqueDist = np.unique(distGrid.flatten())
    ccDist = np.asarray([np.median(ccXY[distGrid == dist]) for dist in uniqueDist])
    idx = np.argsort(ccDist)
    f = PchipInterpolator(ccDist[idx], uniqueDist[idx])
    maxZShift = f(dissim[0, -1]) * 2

    # assign Z positions to each group assuming even spacing across groupZ
    zSpacing = maxZShift / (nZShifts - 1)
    grpZPos = (np.asarray(range(0, nZShifts)) - ((nZShifts - 1) / float(2))) * zSpacing

    newShifts = regDict[shifts]
    zshifts = np.ones(newShifts.shape[0:2])
    for i in range(0, nZShifts):
        for j in grpZIdx[i]:
            zshifts[finalGroups == (j + 1), :] = grpZPos[i]
    finalShifts = np.dstack((newShifts, zshifts))

    regDict['groupZ'] = groupZ
    regDict['grpZIdx'] = grpZIdx
    regDict['grpZPos'] = grpZPos
    regDict['finalShifts'] = finalShifts

    # display results
    ax1 = plt.subplot(2, 1, 1)
    _, bins, _ = ax1.hist(regDict['finalShifts'][:, 1, 2], bins=nZShifts)
    ax1.set_ylabel('Timepoints', color='b')
    ax1.set_xlabel('Z position')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    binCenters = 0.5 * (bins[1:] + bins[:-1])
    ax2 = ax1.twinx()
    Mean = regDict['groupZ'].mean(axis=(1, 2, 3))
    Std = regDict['groupZ'].std(axis=(1, 2, 3))
    ax2.errorbar(binCenters, Mean, Std, ecolor='r', fmt='ro')

    plt.grid(False)
    ax2.set_ylabel(u'Intensity (±SD)', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    plt.show()
    logger.info('Z positions: %s' % (regDict['grpZPos'],))
