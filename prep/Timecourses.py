# Encoding: utf-8
""" Module to interact with the Matlab GUI for spine segmentation 
and pulling out time courses from the data
"""
from __future__ import print_function

import copy
import logging
import os
import sys
import itertools
from builtins import zip
from builtins import map

import SimpleITK as sitk
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objs as go
import plotly.offline as py
import pandas as pd
import scipy.ndimage as ndimg
import skimage.external.tifffile as tf
import tables
import tqdm
from future.utils import iteritems
from networkx.algorithms.dag import descendants
from IPython.display import clear_output
from neurom.io import load_data
from scipy.interpolate import griddata
from scipy.ndimage.morphology import binary_dilation
from scipy.spatial import KDTree
from sklearn import linear_model

from prep.IO import loadmat, writeTiff
from prep.Log import add_logging_to_file
from prep.Utils import getStructure, log_progress, convert_8bit, angle

# Setup logging
logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def initTimeDict(session, maskName='mask', um3bin=100, dilate1=0.5, dilate2=0.5, zStep=1, anatomyZstep=1.6):
    """ Initialize a time courses dictionary

    :param session: current session object
    :param maskName: mask name
    :param um3bin: microns^3 per bin in dendrite
    :param dilate1: first dilation of each spine to subtract from dendrite mask
    :param dilate2: second dilation of the combined mask
    :param zStep: zStep of extended stack
    :param anatomyZstep: Anatomical stack of FOV z spacing
    :return: a timeDict dictionary
    """

    timeDict = dict()
    timeDict['path'] = session.path
    timeDict['fullPath'] = session.path + maskName
    timeDict['fieldsTform'] = session.embedDict['TFormExp']['fieldsTform']
    timeDict['pixelSize'] = np.array([session.pixSizeXY / 1000, session.pixSizeXY / 1000, zStep])
    timeDict['um3bin'] = um3bin
    timeDict['dilate1'] = dilate1
    timeDict['dilate2'] = dilate2
    timeDict['anatomyZstep'] = anatomyZstep
    timeDict['Fs'] = session.volRate
    timeDict['xyPixNum'] = session.Sp['OptionsStruct']['RefPixels']
    timeDict['UM_1X'] = float(session.Sp['OptionsStruct']['UM_1X'])
    timeDict['fieldMask'] = session.fieldMaskBool
    timeDict['imagingLines'] = session.ySize
    timeDict['imagingPixels'] = session.xSize
    timeDict['flyLines'] = session.optFlyLines
    timeDict['nExp'] = session.regDict['nZShifts']
    timeDict['finalShifts'] = session.regDict['finalShifts']
    timeDict['grpZPos'] = session.regDict['grpZPos']
    timeDict['groupZ'] = session.regDict['groupZ']
    timeDict['fix_spines'] = False
    # path to database
    FOV = 'FOV' + str(int(session.meta['FOV'][0]))
    path = os.path.join(session.basePath, 'Database', session.animalID, FOV, session.date + session.run)
    timeDict['databasePath'] = path
    timeDict['SessionPath'] = os.path.join(path, 'Session.tif')
    # add_logging_to_file('prep.Timecourses', timeDict['fullPath'], 'Timecourses.log')
    logger.info('timeDict initialized with full path:' + timeDict['fullPath'] + '.mat')

    return timeDict


def load_hdf5(filename):
    """ loads new .mat file save using the -v7.3 flag in matlab.
    
    :param filename: file to load
    :return: dict with the values
    """
    fh = tables.open_file(filename)
    data = dict()
    data['labelimg'] = fh.root.newCell.labelimg[:].T
    data['dendTable'] = fh.root.newCell.dendTable[:].T
    dend = fh.root.newCell.dend[:]
    dendFix = []
    for d in dend:
        dendFix.append(d[0].T)
    data['dend'] = dendFix
    data['quality'] = fh.root.newCell.quality[:]
    data['note'] = fh.root.newCell.note[:]
    data['dendNum'] = np.squeeze(fh.root.newCell.dendNum[:]).astype(np.uint8)
    return data


def loadMask(timeDict):
    """ load a mask from a mat file

    :param timeDict: time course dictionary needs fullPath
    :return adds: all the data from the mat file
    """
    try:
        logger.info('Trying to load: %s' % timeDict['fullPath'])
        data = loadmat(timeDict['fullPath'][:-4])
    except:
        try:
            data = loadmat(timeDict['fullPath'])
        except (NotImplementedError, TypeError, ValueError):
            logger.info('Trying to load hdf5 at: %s' % timeDict['fullPath'])
            # new -v7.3 mat file is a hdf5 file
            data = load_hdf5(timeDict['fullPath'])

    timeDict['labelimg'] = data['labelimg']
    timeDict['dendTable'] = data['dendTable']
    a = load_data(os.path.join(timeDict['databasePath'], 'dendrite.swc'))
    timeDict['dendR'] = a.data_block[:, 3]
    timeDict['dend'] = data['dend']
    timeDict['dendNumAll'] = data['dendNum'] - 1
    timeDict['quality'] = data['quality']
    timeDict['notes'] = data['note']
    timeDict['dims'] = timeDict['labelimg'].shape
    timeDict['spineMasks'] = len(np.unique(timeDict['labelimg'])) - 1
    timeDict['dendNum'] = len(timeDict['dend'])
    logger.info('imported mask with dims: ' + str(timeDict['labelimg'].shape) + ' ,branches: ' +
                str(timeDict['dendNum']) + ', spine masks:' + str(timeDict['spineMasks']))


def binDendrite(sc, timeDict):
    """ bin the dendrite in um2bin steps

    :param sc: Spark Context
    :param timeDict: time course dictionary needs data from mat file and um2bin param
    :return: adds: dendLabelImg, dendLabelTable
    """
    dims = timeDict['dims']
    labelimg = timeDict['labelimg']
    pixelSize = timeDict['pixelSize']
    dend = timeDict['dend']
    dendTable = timeDict['dendTable']
    dendLabelImg = np.zeros(dims)
    dendLabelTable = np.array([], ndmin=2)
    # last dendrite mask
    counter = np.max(labelimg) + 1
    dendLabelImg = dendLabelImg.reshape(-1)
    dendPathLength = []
    flagFirst = True
    # for each dendrite

    segmentSize = list()
    for i in log_progress(range(timeDict['dendNum']), name='Dendrites'):
        if dend[i].shape[0] == 0:
            logger.info('Skipped dendrite %d' % (i + 1,))
            continue
        X = dend[i][:, 1] - 1
        Y = dend[i][:, 0] - 1
        Z = dend[i][:, 2] - 1
        # find radii
        currentDend = dendTable[dendTable[:, 7] == i + 1, :]
        if len(dend[i].shape) < 2:
            r = currentDend[:, 5] / 2
        else:
            r = griddata(currentDend[:, [2, 3, 4]] / 2, currentDend[:, 5] / 2, dend[i][:, [0, 1, 2]], 'nearest')
        # path length distance
        distPath = np.sqrt((np.diff(X.astype(np.float32) * pixelSize[0])) ** 2 +
                           (np.diff(Y.astype(np.float32)) * pixelSize[1]) ** 2 +
                           (np.diff(Z.astype(np.float32)) * pixelSize[2]) ** 2)
        # find the distance along the dendrite and segment into um2bin parts (10um default)
        dist = np.sqrt((np.diff(X.astype(np.float32))) ** 2 +
                       (np.diff(Y.astype(np.float32))) ** 2 +
                       (np.diff(Z.astype(np.float32))) ** 2)
        dist = dist * ((4.0 / 3.0) * np.pi * r[:-1])
        stopPt = 0
        # logger.info('Dendrite: %d' % (i + 1))
        sys.stdout.flush()
        # for each dendrite segment dilate and add to label image
        # for each dendrite segment dilate and add to label image
        if len(X) == 2:
            starts = [0]
            ends = [2]
            lengths = [np.sum(distPath[0:2])]
        else:
            starts = []
            ends = []
            lengths = []
            while stopPt < len(X):
                startPt = stopPt
                stop1 = np.where(np.cumsum(dist[startPt:]) > timeDict['um3bin'])[0] + startPt
                if len(stop1) > 0:
                    stopPt = np.min(np.array([stop1[0], len(X)]))
                else:
                    stopPt = len(X)
                if startPt == stopPt:
                    stopPt += 2
                starts.append(startPt)
                ends.append(stopPt)
                d = np.sum(distPath[startPt:stopPt])
                if d == 0.0:
                    logger.error('Distance is 0!')
                    break
                else:
                    lengths.append(d)
        rdd = sc.parallelize(list(zip(range(len(starts)), starts, ends)))

        def mask(kse):
            key, start, stop = kse
            dendLabelImgTemp = np.zeros(dims)
            x1 = X[start:stop].astype(int)
            y1 = Y[start:stop].astype(int)
            z1 = Z[start:stop].astype(int)
            # make sure no overflow
            x1[x1 >= dims[0]] = dims[0] - 1
            y1[y1 >= dims[1]] = dims[1] - 1
            z1[z1 >= dims[2]] = dims[2] - 1
            dendLabelImgTemp[x1, y1, z1] = 1
            sSize = np.mean(r[start:stop])
            dendLabelImgTemp = binary_dilation(dendLabelImgTemp > 0,
                                               getStructure(np.array([pixelSize[0], pixelSize[1], 2.5]), sSize))
            index_inner = np.ravel_multi_index(np.where(dendLabelImgTemp > 0), dims)
            return key, sSize, index_inner

        indexRdd = rdd.map(mask).collect()
        for segNum, size, index in indexRdd:
            dendLabelImg[index] = counter
            if flagFirst:
                dendLabelTable = np.array([counter, i, segNum], ndmin=2)
                flagFirst = False
            else:
                dendLabelTable = np.append(dendLabelTable, np.array([counter, i, segNum], ndmin=2), axis=0)
            counter += 1
            segmentSize.append(size)
        dendPathLength.append(np.array(lengths))

    dendLabelImg = dendLabelImg.reshape(labelimg.shape)
    timeDict['dendLabelImg'] = dendLabelImg
    timeDict['dendLabelTable'] = dendLabelTable
    timeDict['dendPathLength'] = dendPathLength
    timeDict['dendSize'] = segmentSize
    writeTiff(timeDict['path'], dendLabelImg, 'dendLabelImg')


def binDendriteByDistance(sc, timeDict, um4bin, hasSoma=True):
    """ bin the dendrite in um4bin steps

    :param sc: Spark Context
    :param timeDict: time course dictionary needs data from mat file and um2bin param
    :return: adds: dendLabelImg, dendLabelTable
    """
    dims = timeDict['dims']
    labelimg = timeDict['labelimg']
    pixelSize = timeDict['pixelSize']
    dend = timeDict['dend']
    dendTable = timeDict['dendTable']
    dendLabelImg = np.zeros(dims)
    dendLabelTable = np.array([], ndmin=2)
    # last dendrite mask
    counter = np.max(labelimg) + 1
    dendLabelImg = dendLabelImg.reshape(-1)
    dendPathLength = []
    flagFirst = True
    # for each dendrite

    if hasSoma:
        if 'cellIndexAll' in timeDict:
            logger.info('Using cellIndexAll')
            soma_dends = np.unique(timeDict['Masks'].loc[timeDict['cellIndexAll']].DendNum.values)
        else:
            logger.info('Using cellIndex')
            soma_dends = np.unique(timeDict['Masks'].loc[timeDict['cellIndex']].DendNum.values)
        logger.info('Found %s as soma dendrite number(s)' % soma_dends)
    else:
        soma_dends = []
    valid = [len(dend[k].shape) == 2 for k in soma_dends]
    soma_dends = soma_dends[valid]
    timeDict['soma_dends'] = soma_dends
    segmentSize = list()
    for i in log_progress(range(timeDict['dendNum']), name='Dendrites'):
        if (dend[i].shape[0] == 0) | (len(dend[i].shape) < 2):
            logger.info('Skipped dendrite %d' % (i + 1,))
            continue
        if len(soma_dends)>0:
            if i == soma_dends[0]:
                somaFlag = True
            else:
                somaFlag = False
        else:
            somaFlag = False
        if len(soma_dends) >= 2 and somaFlag:
            # soma mask needs stiching
            a = [(dend[k][:, 1] - 1) for k in soma_dends]
            X = np.hstack(a)
            a = [(dend[k][:, 0] - 1) for k in soma_dends]
            Y = np.hstack(a)
            a = [(dend[k][:, 2] - 1) for k in soma_dends]
            Z = np.hstack(a)
            indexs = [dendTable[:, 7] == k + 1 for k in soma_dends]
            index = np.logical_or(indexs[0], indexs[1])
            currentDend = dendTable[index, :]
        else:
            if len(soma_dends) >= 2 and i in soma_dends:
                continue
            X = dend[i][:, 1] - 1
            Y = dend[i][:, 0] - 1
            Z = dend[i][:, 2] - 1
            # find radii
            currentDend = dendTable[dendTable[:, 7] == i + 1, :]
        if len(dend[i].shape) < 2:
            r = currentDend[:, 5] / 2
        else:
            r = griddata(currentDend[:, [2, 3, 4]] / 2, currentDend[:, 5] / 2, dend[i][:, [0, 1, 2]], 'nearest')
        # path length distance
        dist = np.sqrt((np.diff(X.astype(np.float32) * pixelSize[0])) ** 2 +
                       (np.diff(Y.astype(np.float32)) * pixelSize[1]) ** 2 +
                       (np.diff(Z.astype(np.float32)) * pixelSize[2]) ** 2)

        if somaFlag:
            starts = [0]
            ends = [len(X)]
            lengths = [np.sum(dist)]
        else:
            if len(X) == 2:
                pass  # continue
            else:
                stopPt = 0
                starts = []
                ends = []
                lengths = []
                while stopPt < len(X):
                    startPt = stopPt
                    stop1 = np.where(np.cumsum(dist[startPt:]) > um4bin)[0] + startPt
                    if len(stop1) > 0:
                        stopPt = np.min(np.array([stop1[0], len(X)]))
                    else:
                        break
                    starts.append(startPt)
                    ends.append(stopPt)
                    d = np.sum(dist[startPt:stopPt])
                    if d == 0.0:
                        logger.error('Distance is 0!')
                        break
                    else:
                        lengths.append(d)
        rdd = sc.parallelize(list(zip(range(len(starts)), starts, ends)))

        def mask(kse):
            key, start, stop = kse
            dendLabelImgTemp = np.zeros(dims)
            x1 = X[start:stop].astype(int)
            y1 = Y[start:stop].astype(int)
            z1 = Z[start:stop].astype(int)
            # make sure no overflow
            x1[x1 >= dims[0]] = dims[0] - 1
            y1[y1 >= dims[1]] = dims[1] - 1
            z1[z1 >= dims[2]] = dims[2] - 1
            dendLabelImgTemp[x1, y1, z1] = 1
            sSize = np.mean(r[start:stop])
            dendLabelImgTemp = binary_dilation(dendLabelImgTemp > 0,
                                               getStructure(np.array([pixelSize[0], pixelSize[1], 2.5]), sSize))
            index_inner = np.ravel_multi_index(np.where(dendLabelImgTemp > 0), dims)
            return key, sSize, index_inner

        indexRdd = rdd.map(mask).collect()
        for segNum, size, index in indexRdd:
            if len(index)>0:
                dendLabelImg[index] = counter
                if flagFirst:
                    dendLabelTable = np.array([counter, i, segNum], ndmin=2)
                    flagFirst = False
                else:
                    dendLabelTable = np.append(dendLabelTable, np.array([counter, i, segNum], ndmin=2), axis=0)
                counter += 1
                segmentSize.append(size)
        dendPathLength.append(np.array(lengths))

    dendLabelImg = dendLabelImg.reshape(labelimg.shape)
    timeDict['dendLabelImgB'] = dendLabelImg
    timeDict['dendLabelTableB'] = dendLabelTable
    timeDict['dendPathLengthB'] = dendPathLength
    timeDict['dendSizeB'] = segmentSize
    writeTiff(timeDict['path'], dendLabelImg, 'dendLabelImgB')


def getMasks(timeDict, makeB=False):
    """ get spine masks by dilating twice

    :param timeDict: time course dictionary need dendrite masks and pixelSize, dilate1 and dilate2
    :return: adds labelimgAll
    """
    labelimg = timeDict['labelimg']
    if not makeB:
        dendLabelImg = timeDict['dendLabelImg']
    else:
        dendLabelImg = timeDict['dendLabelImgB']
    pixelSize = timeDict['pixelSize']
    labelimg2 = copy.deepcopy(labelimg)
    labelimg2 = labelimg2.reshape(-1)
    rangeEnd = np.max(labelimg) + 1
    # logger.info('Going over spines')
    for i in log_progress(range(1, rangeEnd), name='Spines'):
        index = np.ravel_multi_index(np.where(labelimg == i), labelimg.shape)
        labelimgTemp = np.zeros(labelimg.shape)
        labelimgTemp = labelimgTemp.reshape(-1)
        labelimgTemp[index] = True
        labelimgTemp = labelimgTemp.reshape(labelimg.shape)
        labelimgTemp = binary_dilation(labelimgTemp, getStructure(pixelSize, timeDict['dilate1']))
        index2 = np.ravel_multi_index(np.where(labelimgTemp), labelimg.shape)
        labelimg2[index2] = i
    labelimg2 = labelimg2.reshape(labelimg.shape)
    labelimg3 = np.invert(binary_dilation(labelimg2 > 0, getStructure(pixelSize, timeDict['dilate2'])))
    dendLabelImg2 = dendLabelImg * labelimg3
    labelimgAll = dendLabelImg2 + labelimg
    if not makeB:
        timeDict['labelimgAll'] = labelimgAll
    else:
        timeDict['labelimgAllB'] = labelimgAll


def getMasksDF(timeDict, makeB=False, plot=True):
    """ get Masks DataFrame and fix dendrite assignments

    :param timeDict: time course dictionary
    :param fix: if True will try to reassign spine numbers to dendrite numbers based on distance
    :return: adds the Masks DataFrame
    """
    fix = timeDict['fix_spines']
    # todo: check for nans
    if not makeB:
        labelimgAll = timeDict['labelimgAll']
        dendLabelTable = timeDict['dendLabelTable']
        dendPathLength = timeDict['dendPathLength']
    else:
        labelimgAll = timeDict['labelimgAllB']
        dendLabelTable = timeDict['dendLabelTableB']
        dendPathLength = timeDict['dendPathLengthB']

    mask = (labelimgAll > 0).astype(int)
    label = labelimgAll
    nMasks = label.max().astype(int)
    centers = np.array(ndimg.center_of_mass(mask, label, range(1, nMasks + 1)))
    pixelSize = timeDict['pixelSize'][0]
    centers[:, 0:2] = centers[:, 0:2] * pixelSize
    Centers = pd.Series([tuple(row) for row in centers])
    types = ['Spine'] * timeDict['dendNumAll'].shape[0] + ['Dendrite'] * dendLabelTable.shape[0]
    MaskType = pd.Series(types)
    index = dendLabelTable[:, 0]
    DendNum = pd.Series(np.concatenate((timeDict['dendNumAll'], dendLabelTable[:, 1].T)))
    DendSegment = pd.Series(dendLabelTable[:, 2], index=index - 1)
    Masks = pd.DataFrame({'Centers': Centers,
                          'MaskType': MaskType,
                          'DendNum': DendNum,
                          'DendSegment': DendSegment})
    dendrites = Masks[(Masks['MaskType'] == 'Dendrite')]['DendNum'].unique()
    dendCenters = dict()
    for dend in dendrites:
        dendCenters[dend] = np.asarray(
            list(Masks[(Masks['MaskType'] == 'Dendrite') & (Masks['DendNum'] == dend)]['Centers'].values))
    spines = Masks[(Masks['MaskType'] == 'Spine')]['DendNum'].unique()
    if fix:
        indexReplace = [None] * len(spines)
        for index, S_dend in enumerate(spines):
            indexReplace[index] = Masks[(Masks['MaskType'] == 'Spine') & (Masks['DendNum'] == S_dend)].index
        new_values = [None] * len(spines)
        for index, S_dend in enumerate(spines):
            SpineCenters = np.asarray(
                list(Masks[(Masks['MaskType'] == 'Spine') & (Masks['DendNum'] == S_dend)]['Centers'].values))
            minDist = np.zeros((SpineCenters.shape[0], dendrites.shape[0]))
            for indexS, spine in enumerate(SpineCenters):
                for indexD, dend in enumerate(dendrites):
                    currentDend = dendCenters[dend]
                    minDist[indexS, indexD] = np.min(np.sum((currentDend - spine) ** 2, axis=1))
            # indexReplace = Masks[(Masks['MaskType'] == 'Spine') & (Masks['DendNum'] == S_dend)].index
            new_values[index] = dendrites[np.argmin(minDist.mean(axis=0))].astype(int)
            # Masks['DendNum'][indexReplace[index]] = new_values
            logger.info('Old: %d, new: %d' % (S_dend, new_values[index]))
        for index, S_dend in enumerate(spines):
            Masks['DendNum'][indexReplace[index]] = new_values[index]
    Masks.dendPathLength = np.hstack(dendPathLength)
    if not makeB:
        timeDict['Masks'] = Masks
    else:
        for index in Masks[(Masks['MaskType'] == 'Spine')].index:
            Masks['DendNum'][index] = timeDict['Masks']['DendNum'][index]
        timeDict['MasksB'] = Masks

    if plot:
        SpineCenters = np.asarray(list(Masks[(Masks.MaskType == 'Spine')]['Centers']))
        SpineDend = np.asarray(list(Masks[(Masks.MaskType == 'Spine')]['DendNum']))
        plt.figure(figsize=(13, 13))
        plt.imshow(labelimgAll.max(axis=2).transpose(1, 0))
        for center, dend in zip(SpineCenters, SpineDend):
            plt.text(center[0] / pixelSize, center[1] / pixelSize, str(dend))
        dendNums = np.unique(dendLabelTable[:, 1].T)
        for num in dendNums:
            current = Masks[(Masks.MaskType == 'Dendrite') & (Masks.DendNum == num)].Centers
            XY = (np.array(list(map(list, current.values)))[:, [0, 1]].mean(axis=0) / pixelSize).astype(int)
            plt.text(XY[0], XY[1], str(num), size=20, color='white')


def getFieldMasks(timeDict, makeB=False):
    """ assigns the masks back to the original data format

    :param timeDict: time course dictionary need labelimgAll, fieldsTform
    :return: adds: labelList
    """
    if not makeB:
        labelimgAll = timeDict['labelimgAll']
    else:
        labelimgAll = timeDict['labelimgAllB']
    fieldsTform = timeDict['fieldsTform']
    linearLabelimg = labelimgAll.flatten(order='F')
    nLabel = np.max(linearLabelimg).astype(np.int64)
    labelInfoAll = [None] * len(fieldsTform)
    for i in range(0, len(fieldsTform)):
        Tform = fieldsTform[i].T
        labelInfo = np.zeros((Tform.shape[0], 3))
        for j in range(0, Tform.shape[0]):
            try:
                labelInfo[j, 0] = linearLabelimg[int(Tform[j, 0])]
            except Exception:
                logger.error(str(i) + ' ' + str(j))
                labelInfo[j, 0] = linearLabelimg[int(Tform[j, 0])]
            labelInfo[j, 1] = Tform[j, 1]
            labelInfo[j, 2] = i
        labelInfoAll[i] = labelInfo.T

    labelArray = np.hstack(labelInfoAll)
    labelList = [None] * nLabel
    for i in range(0, nLabel):
        labelList[i] = labelArray[1:3, labelArray[0, :] == (i + 1)]

    if not makeB:
        timeDict['labelList'] = labelList
    else:
        timeDict['labelListB'] = labelList


def getRegionData(sc, data, timeDict, makeB=False):
    """ get time course data

    :param sc: SparkContext
    :param data: Thunder Images object
    :param timeDict: time course dictionary needs labelList
    :return: adds: TC, TCPixels, TCMotion
    """
    # todo: if tracing out of bounds!
    if not makeB:
        labelList = copy.deepcopy(timeDict['labelList'])
    else:
        labelList = copy.deepcopy(timeDict['labelListB'])
    # maxIndex = np.prod(data.shape[1:]) * len(grpZPos)
    # for i, label in enumerate(labelList):
    #     index = label[0, :] < maxIndex
    #     labelList[i][0, :] =labelList[i][0, :][index]
    #     index = label[1, :] < maxIndex
    #     labelList[i][1, :] =labelList[i][1, :][index]

    finalShifts = timeDict['finalShifts']
    grpZPos = timeDict['grpZPos']
    groupZ = timeDict['groupZ']
    labelListBC = sc.broadcast(labelList)
    shiftsBC = sc.broadcast(finalShifts)
    RegMean = groupZ.transpose(1, 2, 0, 3).flatten(order='F')
    RegMeanBC = sc.broadcast(RegMean)

    # old = np.seterr(all='raise')
    # logger.info('Set numpy errors to raise')

    def offsetVol(vol, pos, grpZPos2):
        np.seterr(all='raise')
        out = np.zeros((vol.shape[0], vol.shape[1], grpZPos2.shape[0] * vol.shape[2]))
        out[:] = np.NAN
        offset = np.argmin(np.absolute(grpZPos2 - pos))
        for i in range(0, vol.shape[2]):
            newFieldID = i * grpZPos2.shape[0] + offset
            out[:, :, newFieldID] = vol[:, :, i]
        return out.flatten(order='F').astype('float32')

    def offsetVolPar(kv):
        np.seterr(all='raise')
        key, ary = kv
        return offsetVol(ary, shiftsBC.value[np.array(key).astype(int), 0, 2], grpZPos).astype('float32')

    def nanMeanByRegions(kv):
        np.seterr(all='raise')
        key, ary = kv
        mean_values = []
        for grp in labelListBC.value:
            a = np.nansum(ary[np.array(grp[1, :], dtype='uint32')] * grp[0, :].flatten())
            b = np.sum(~np.isnan(ary[np.array(grp[1, :], dtype='uint64')]) * grp[0, :].flatten())
            if b == 0.0:
                mean_values.append(np.nan)
            else:
                mean_values.append(a / b)
        # mean_values = [np.nansum(ary[np.array(grp[1, :], dtype='uint32')] * grp[0, :].flatten()) / np.sum(
        #     ~np.isnan(ary[np.array(grp[1, :], dtype='uint64')]) * grp[0, :].flatten()) for grp in labelListBC.value]
        return np.array(mean_values, dtype=ary.dtype).reshape((1, 1, -1))

    def nanMeanByRegionsMotion(kv):
        np.seterr(all='raise')
        key, ary = kv
        compMean = RegMeanBC.value * np.absolute(np.sign(ary))
        norm_values = []
        for grp in labelListBC.value:
            a = np.nansum(compMean[np.array(grp[1, :], dtype='uint32')] * grp[0, :].flatten())
            b = np.sum(~np.isnan(compMean[np.array(grp[1, :], dtype='uint32')]) * grp[0, :].flatten())
            if b == 0.0:
                norm_values.append(np.nan)
            else:
                norm_values.append(a / b)
        return np.array(norm_values, dtype='float32').reshape((1, 1, -1))

    def nanMeanByRegionsPixels(kv):
        np.seterr(all='raise')
        key, ary = kv
        pixels = [np.sum(~np.isnan(ary[np.array(grp[1, :], dtype='uint64')]) * grp[0, :].flatten()) for grp in
                  labelListBC.value]
        return np.array(pixels, dtype='float32').reshape((1, 1, -1))

    RegDataExp = data.map(offsetVolPar, with_keys=True)
    RegDataExp.cache()
    RegDataExp.count()
    if not makeB:
        timeDict['TC'] = RegDataExp.map(nanMeanByRegions, with_keys=True).toarray().T
        logger.info('Got TC')
        timeDict['TCMotion'] = RegDataExp.map(nanMeanByRegionsMotion, with_keys=True).toarray().T
        logger.info('Got TCMotion')
        timeDict['TCPixels'] = RegDataExp.map(nanMeanByRegionsPixels, with_keys=True).toarray().T
        logger.info('Got TCPixels')
    else:
        timeDict['TCB'] = RegDataExp.map(nanMeanByRegions, with_keys=True).toarray().T
        logger.info('Got TCB')
        timeDict['TCMotionB'] = RegDataExp.map(nanMeanByRegionsMotion, with_keys=True).toarray().T
        logger.info('Got TCMotionB')
        timeDict['TCPixelsB'] = RegDataExp.map(nanMeanByRegionsPixels, with_keys=True).toarray().T
        logger.info('Got TCPixelsB')
    RegDataExp.uncache()
    # np.seterr(**old)
    # logger.info('Set numpy errors to: %s' % old)


def getBaseline(sc, timeDict, regWindow=1000, maxZscore=2.0, step=8, makeB=False):
    """ estimates baseline

    :param sc: Spark Context
    :param timeDict: time course dictionary needs: TC and TCMotion
    :param regWindow: number of time points to smooth
    :param maxZscore: maximum z score for a good point
    :param step: 
    :return: adds: TCBaseline
    """

    from sklearn import linear_model
    # t = time.time()
    timeDict['regWindow'] = regWindow
    timeDict['maxZscore'] = maxZscore
    timeDict['step'] = step

    def estBaseline(key, model_inner):
        np.seterr(all='warn')
        from scipy.stats import gaussian_kde
        start = int(max((0, key - regWindow / 2)))
        stop = int(min((len(x), key + regWindow / 2)))
        x1 = x_BC.value[start:stop]
        y1 = y_BC.value[start:stop]
        # p1 = p[start:stop]
        x2 = x1[np.logical_not(np.isnan(x1))]
        y2 = y1[np.logical_not(np.isnan(y1))]
        if np.any(y2):
            if len(y2) > 100:
                kernel = gaussian_kde(y2)
                low, high = np.percentile(y2, [25, 75]).astype(int)
                step_inner = (high - low) / 100.
                testRange = low + np.arange(start=1, stop=101, dtype=int) * step_inner
                estMode = testRange[np.argmax(kernel(testRange))]
            else:
                estMode = np.median(y2)

            y3 = y2[(y2 - estMode) < 0] - estMode
            std = np.std(np.hstack((y3, -y3)))
            zscore = (y1 - estMode) / std
            goodPts = np.logical_and((zscore < maxZscore), not_nan_BC.value[start:stop])
        else:
            goodPts = []
        if np.any(goodPts):
            model_inner = model_inner.fit(x1[goodPts].reshape(-1, 1), y1[goodPts].reshape(-1, 1))
            coef_inner = model_inner.coef_
            if coef_inner < 0.1:
                coef_inner = np.nanmean(y2) / np.nanmean(x2)
        else:
            coef_inner = np.NAN
        return key, np.squeeze(coef_inner)

    model = linear_model.LinearRegression(fit_intercept=False)
    if not makeB:
        TC = timeDict['TC']
        TCMotion = timeDict['TCMotion']
        TCPixels = timeDict['TCPixels']
    else:
        TC = timeDict['TCB']
        TCMotion = timeDict['TCMotionB']
        TCPixels = timeDict['TCPixelsB']
    TCBaselineDict = dict()
    inter_x = np.arange(0, TC.shape[1], 1, int)
    inter_xp = np.arange(0, TC.shape[1], step, int)
    for i in tqdm.tqdm(range(TC.shape[0])):
        x = TCMotion[i, :]
        y = TC[i, :]
        p = TCPixels[i, :]
        not_nan = np.logical_and(np.logical_not(np.isnan(x)), np.logical_not(np.isnan(y)))
        not_nan = np.logical_and(not_nan, np.logical_not(np.isnan(p)))
        x_BC = sc.broadcast(x)
        y_BC = sc.broadcast(y)
        not_nan_BC = sc.broadcast(not_nan)
        coefDict = sc.parallelize(range(0, len(x), step)).map(lambda x2: estBaseline(x2, model)).collectAsMap()
        coef = np.array([coefDict[idx] for idx in range(0, len(x), step)])
        coef = np.interp(inter_x, inter_xp, coef)
        TCBaselineDict[i] = x * np.squeeze(coef)
        x_BC.unpersist()
        y_BC.unpersist()
        not_nan_BC.unpersist()
        # current = time.time() - t
        # m, s = divmod(current, 60)
        # logger.info('i: %d, %02d:%02d' % (i, m, s))
        # sys.stdout.flush()

    TCBaseline = np.array([np.squeeze(TCBaselineDict[idx]) for idx in TCBaselineDict.keys()])
    if not makeB:
        timeDict['TCBaseline'] = TCBaseline
        old = np.seterr(all='warn')
        timeDict['TCdiv'] = TC / TCBaseline - 1
    else:
        timeDict['TCBaselineB'] = TCBaseline
        old = np.seterr(all='warn')
        timeDict['TCdivB'] = TC / TCBaseline - 1


def getNoise(sc, timeDict, makeB=False):
    """ estimate noise

    :param sc: Spark Context
    :param timeDict: time course dictionary needs: TC and TCMotion
    :return: adds TCNoise and TCZscore
    """

    def fitNoise(key):
        from scipy.optimize import curve_fit

        def model_noise(x2, Ndiv, Nscale, offset):
            lambdaAct = ((x2 + offset) ** Nscale) / Ndiv
            return (lambdaAct ** 0.5) / lambdaAct

        idx_inner = (TCdiv_BC.value[key[0], :] < 0).nonzero()[0]
        TCPixels2 = copy.deepcopy(TCPixels_BC.value[key[0], :])
        TCPixels2[TCPixels_BC.value[key[0], :] == 0] = np.NAN
        x = TCPixels2[idx_inner] * TCBaseline_BC.value[key[0], idx_inner]
        y = -TCdiv_BC.value[key[0], idx_inner]
        validTps = np.isfinite(x) & np.isfinite(y)
        x = x[validTps]
        y = y[validTps]
        x2 = x[x > 0]
        y2 = y[x > 0]
        try:
            opt_parameters, parm_cov = curve_fit(model_noise, x2, y2 ** 2, maxfev=10000, method='trf')
            if np.any(np.logical_not(np.isfinite(opt_parameters))):
                TC_noise = np.ones_like(TCPixels2) * np.mean(y2 ** 2) ** 0.5
            else:
                TC_noise = model_noise(TCPixels2 * TCBaseline_BC.value[key[0], :], opt_parameters[0],
                                       opt_parameters[1], opt_parameters[2]) ** 0.5
        except:
            TC_noise = np.ones_like(TCPixels2) * np.mean(y2 ** 2) ** 0.5
        return key, TC_noise

    if not makeB:
        TCdiv_BC = sc.broadcast(timeDict['TCdiv'])
        TCPixels_BC = sc.broadcast(timeDict['TCPixels'])
        TCBaseline_BC = sc.broadcast(timeDict['TCBaseline'])
        idxList = list(zip(range(0, timeDict['TCdiv'].shape[0])))
    else:
        TCdiv_BC = sc.broadcast(timeDict['TCdivB'])
        TCPixels_BC = sc.broadcast(timeDict['TCPixelsB'])
        TCBaseline_BC = sc.broadcast(timeDict['TCBaselineB'])
        idxList = list(zip(range(0, timeDict['TCdivB'].shape[0])))

    fitNoiseDict = sc.parallelize(idxList).map(fitNoise).collectAsMap()
    TCNoise = np.array([fitNoiseDict[idx] for idx in idxList])

    if not makeB:
        timeDict['TCNoise'] = TCNoise
        timeDict['TCzscore'] = timeDict['TCdiv'] / TCNoise
    else:
        timeDict['TCNoiseB'] = TCNoise
        timeDict['TCzscoreB'] = timeDict['TCdivB'] / TCNoise


def loadInfo(timeDict):
    """ load the .mat file with session info

    :param timeDict: time course dictionary
    :return: loaded mat file
    """
    Info = loadmat(os.path.join(timeDict['databasePath'], 'prepareMasksAutoSWC.mat'))
    timeDict['InfoSWC'] = Info['Info']
    logger.info('Loaded ' + timeDict['path'] + 'prepareMasksAutoSWC.mat')
    Info = loadmat(os.path.join(timeDict['databasePath'], 'prepareMasksAuto.mat'))
    timeDict['Info'] = Info['Info']
    logger.info('Loaded ' + timeDict['path'] + 'prepareMasksAuto.mat')


def getTransform(timeDict, outputFile='InvTransform.h5', getAligned=True, do_8bit=True, sat=1):
    """ calculates a affine transformation from session space to anatomy stack space

    :param timeDict: time course dict
    :param outputFile: name of transform file to save to
    :param getAligned: flag to apply the transformation and return 'Aligned' to timeDict
    :return: inverse transformation (from session space to anatomical stack space)
    """

    # callback invoked when the StartEvent happens, sets up our new data
    def start_plot():
        global metric_values, multires_iterations

        metric_values = []
        multires_iterations = []

    # callback invoked when the EndEvent happens, do cleanup of data and figure
    def end_plot():
        global metric_values, multires_iterations
        try:
            del metric_values
            del multires_iterations
            # close figure, we don't want to get a duplicate of the plot latter on
            plt.close()
        except Exception:
            pass

    # callback invoked when the IterationEvent happens, update our data and display new figure
    def plot_values(registration_method_inner):
        global metric_values, multires_iterations
        val = registration_method_inner.GetMetricValue()
        if np.isfinite(val):
            metric_values.append(registration_method_inner.GetMetricValue())
            # clear the output area (wait=True, to reduce flickering), and plot current data
            clear_output(wait=True)
            # plot the similarity metric values
            plt.plot(metric_values, 'r')
            plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
            plt.xlabel('Iteration Number', fontsize=12)
            plt.ylabel('Metric Value', fontsize=12)
            plt.show()

    # callback invoked when the sitkMultiResolutionIterationEvent happens, update the index into the
    # metric_values list.
    def update_multires_iterations():
        global metric_values, multires_iterations
        multires_iterations.append(len(metric_values))

    pxSize = timeDict['pixelSize'][0]
    xyStep = timeDict['UM_1X'] / timeDict['xyPixNum']
    aMaskSmooth = tf.imread(timeDict['SessionPath']).transpose(1, 2, 0)
    sMaskSmooth = tf.imread(os.path.join(timeDict['databasePath'], 'expended_new.tif')).transpose(2, 1, 0)
    sMaskSmooth2 = sMaskSmooth.flatten()[np.logical_not(np.isnan(sMaskSmooth.flatten()))]
    t = np.median(sMaskSmooth2)
    sMaskSmooth[sMaskSmooth < t] = 0
    sMaskSmooth = np.nan_to_num(sMaskSmooth)
    zNum = aMaskSmooth.shape[2]
    timeDict['zNum'] = zNum
    timeDict['sMaskSmooth'] = sMaskSmooth
    timeDict['aMaskSmooth'] = aMaskSmooth
    logger.info('Prepared images')
    sys.stdout.flush()

    target = sitk.GetImageFromArray(aMaskSmooth)
    moving = sitk.GetImageFromArray(sMaskSmooth)
    target.SetSpacing((timeDict['anatomyZstep'], xyStep, xyStep))
    moving.SetSpacing((timeDict['pixelSize'][2], pxSize, pxSize))

    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsCorrelation()
    registration_method.SetInterpolator(sitk.sitkLinear)

    # optimizer settings
    registration_method.SetOptimizerAsGradientDescent(learningRate=0.5, numberOfIterations=20,
                                                      convergenceMinimumValue=1e-6, convergenceWindowSize=10,
                                                      estimateLearningRate=sitk.ImageRegistrationMethod.EachIteration)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # setup for the multi-resolution framework
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # connect all of the observers so that we can perform plotting during registration
    registration_method.AddCommand(sitk.sitkStartEvent, start_plot)
    registration_method.AddCommand(sitk.sitkEndEvent, end_plot)
    registration_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, update_multires_iterations)
    registration_method.AddCommand(sitk.sitkIterationEvent, lambda: plot_values(registration_method))

    transform = sitk.CenteredTransformInitializer(target, moving, sitk.AffineTransform(3),
                                                  sitk.CenteredTransformInitializerFilter.MOMENTS)
    registration_method.SetInitialTransform(transform)

    registration_method.Execute(target, moving)
    logger.info('Final metric value: {0}'.format(registration_method.GetMetricValue()))
    logger.info('Optimizer\'s stopping condition, {0}'.format(
        registration_method.GetOptimizerStopConditionDescription()))
    if getAligned:
        reSampler = sitk.ResampleImageFilter()
        reSampler.SetTransform(transform)
        reSampler.SetSize(target.GetSize())
        reSampler.SetOutputSpacing(target.GetSpacing())
        aligned = sitk.GetArrayFromImage(reSampler.Execute(moving))
        timeDict['aligned'] = aligned
        path = os.path.join(timeDict['path']+'View', '')
        writeTiff(path, aligned, 'aligned')
        if do_8bit:
            aligned_8bit = convert_8bit(aligned, sat_percent=sat, ignore_nans=False, ignore_zero=True)
            writeTiff(path, aligned_8bit, 'aligned_8bit', dtype='uint8')


    transformInv = transform.GetInverse()
    sitk.WriteTransform(transformInv, timeDict['path'] + outputFile)
    logger.info('Saved inverse transform to ' + timeDict['path'] + outputFile)
    return transformInv


def transformPoints(timeDict, transformInv, useSWC=True, makeB=False):
    """ transform the center points of all masks to anatomy space in um

    :param timeDict: time course dictionary
    :param transformInv: the SimpleITK transformation
    :param useSWC: load data from swc reconstruction file
    :return: adds AnatomyCenters, ConnectingPoint and ConnectingDist to Masks DataFrame
    """
    # transform point to anatomy space
    if useSWC:
        Info = timeDict['InfoSWC']
        All = Info['AllSWC']
    else:
        Info = timeDict['Info']
        All = Info['All']
    img = tf.imread(timeDict['SessionPath']).astype(int).transpose(1, 2, 0)
    imgSize = img.shape
    xyStep = timeDict['UM_1X'] / timeDict['xyPixNum']
    if not makeB:
        Masks = timeDict['Masks']
    else:
        Masks = timeDict['MasksB']
    AnatomyCenters = list([])
    for center in Masks.Centers:
        point = transformInv.TransformPoint((center[2], center[0], center[1]))
        AnatomyCenters.append((point[1], point[2], point[0]))
    Masks['AnatomyCenters'] = AnatomyCenters

    # find closest point in interpolated space
    Connecting = list([])
    segList = list([])
    for segNum, segment in enumerate(All):
        x, y, z = np.unravel_index(segment, imgSize, 'F')
        # subtract one for matlab 1 indexed problem
        x = x - 1
        y = y - 1
        z = z - 1
        # move to real space coordinates (um)
        x = x * xyStep
        y = y * xyStep
        z = z * timeDict['anatomyZstep']
        if len(x.shape) > 0:
            Connecting.extend(list(zip(x, y, z)))
            segList.extend([segNum] * x.shape[0])
        else:
            Connecting.append((x, y, z))
            segList.append(segNum)
    Connecting = np.array(Connecting)
    tree = KDTree(Connecting)
    distances, indexes = tree.query(AnatomyCenters)
    loc = Connecting[indexes, :]
    seg = np.array(segList)[indexes]
    loc2 = list(map(tuple, loc))
    Masks['Segment'] = seg
    Masks['ConnectingPoint'] = loc2
    Masks['ConnectingDist'] = distances

    # find closest segment
    if useSWC:
        Table2 = np.array(timeDict['InfoSWC']['TableSWC'])
    else:
        Table2 = np.array(timeDict['Info']['Table'])
    Table = Table2[:, 2:5]
    Table[:, 0] = Table[:, 0]  # * xyStep
    Table[:, 1] = Table[:, 1]  # * xyStep
    Table[:, 2] = Table[:, 2]  # * timeDict['anatomyZstep']
    tree2 = KDTree(Table)
    distances2, indexes2 = tree2.query(AnatomyCenters)
    Ids = Table2[indexes2, 0]
    Masks['ParentId'] = Ids
    # timeDict['Masks'] = Masks


def loadTransform(path, filename='InvTransform.h5'):
    """ loads a SimpleITK transformation object

    :param path: path to file (session.path)
    :param filename: transformation filename (*.h5)
    :return: the transformation
    """

    return sitk.ReadTransform(path + filename)


def getExcludeIndexB(timeDict):
    li=timeDict['labelimgAll']
    liB=timeDict['labelimgAllB']
    eI = timeDict['excludeIndex']
    Masks = timeDict['Masks']
    spineIdx = np.asarray(list(Masks[(Masks.MaskType == 'Spine')].index))
    dendIdx = np.asarray(list(Masks[(Masks.MaskType == 'Dendrite')].index))
    eSpine = np.intersect1d(spineIdx,eI)
    eDend = np.intersect1d(dendIdx,eI)
    if np.any(np.isfinite(eDend)) and len(eDend)>0:
        liMask=np.zeros(li.shape)
        for idx in eDend:
            liMask[li==(idx+1)]=1
        eDendB = np.setdiff1d(np.unique(liB[liMask.astype('bool')]),0)-1
        eDendNum = timeDict['Masks'].DendNum[eDend]
        eDendNumB = timeDict['MasksB'].DendNum[eDendB]
        keepIdx = np.intersect1d(eDendNumB,eDendNum)
        eDendB=eDendB[[np.any(keepIdx==x) for x in eDendNumB]]
        timeDict['excludeIndexB'] = np.hstack([eSpine,eDendB])
    else:
        timeDict['excludeIndexB'] = eSpine
    timeDict['excludeIndexB'] = timeDict['excludeIndexB'].astype(int)


def get_graph_from_swc(timeDict, session):
    """

    :param timeDict:
    :param session:
    :return:
    """
    name = session.Sp['OptionsStruct']['filename']
    name = str(name[name.rfind('\\') + 1:])
    FOV = timeDict['databasePath'][:timeDict['databasePath'].rfind('/') + 1]
    table = load_data(os.path.join(FOV, name))
    logger.info('loaded SWC from %s: ' % os.path.join(FOV, name))
    table2 = pd.DataFrame(data=table.data_block, columns=['x', 'y', 'z', 'r', 'type', 'id', 'pID'])
    table2.z = table2.z * 1.6
    table2.x = table2.x * (timeDict['UM_1X'] / timeDict['xyPixNum'])
    table2.y = table2.y * (timeDict['UM_1X'] / timeDict['xyPixNum'])
    table2 = table2[['id', 'type', 'x', 'y', 'z', 'r', 'pID']]
    weight = []
    start = []
    end = []
    for line in table2.iterrows():
        Id = int(line[0]) + 1
        p1 = np.array(line[1][2:5])
        pID = int(line[1][6])
        if Id == 0 or pID == -1:
            weight.append(0)
            start.append(-1)
            end.append(1)
            continue
        p2 = np.array(table2.loc[int(pID - 1)][2:5])
        length = np.linalg.norm(p2 - p1)
        start.append(pID)
        end.append(Id)
        weight.append(length)
    table3 = pd.DataFrame(data=np.array(list(zip(start, end, weight))), columns=['start', 'end', 'weight'])
    G = nx.from_pandas_edgelist(table3, source='start', target='end', edge_attr=['weight'])
    isTree = nx.algorithms.tree.recognition.is_tree(G)
    logger.info('Is tree: %s' % isTree)
    timeDict['graph'] = G


def get_path_length(timeDict, makeB=False, plot=True):
    G = timeDict['graph']
    pathLen = dict(nx.shortest_path_length(G, weight='weight'))
    CellBodyPath = pathLen[1.0]
    timeDict['CellBodyPath'] = CellBodyPath
    if not makeB:
        Masks = timeDict['Masks']
        labelimgAll = timeDict['labelimgAll']
    else:
        Masks = timeDict['MasksB']
        labelimgAll = timeDict['labelimgAllB']
    for i, mask in enumerate(Masks.ParentId):
        if mask != 1.0:
            Masks.set_value(i, 'PathLength', CellBodyPath[mask])
        else:
            Masks.set_value(i, 'PathLength', 0)
    pathLenArray = np.zeros((len(pathLen) - 1, len(pathLen) - 1))
    for key, value in iteritems(pathLen):
        if key != -1.0:
            for key2 in sorted(value):
                if key2 != -1.0:
                    pathLenArray[int(key - 1), int(key2 - 1)] = value[key2]
    timeDict['pathLengthAll'] = pathLenArray
    if not makeB:
        timeDict['pathIndex'] = np.array(np.argsort(Masks.PathLength))
        timeDict['pathLength'] = np.sort(np.array(Masks.PathLength))
    else:
        timeDict['pathIndexB'] = np.array(np.argsort(Masks.PathLength))
        timeDict['pathLengthB'] = np.sort(np.array(Masks.PathLength))

    if plot:
        plt.figure(figsize=(12, 12))
        plt.imshow(labelimgAll.max(axis=2).transpose(1, 0))
        DendCenters = np.asarray(list(Masks[(Masks.MaskType == 'Dendrite')]['Centers']))
        DendPath = np.asarray(list(Masks[(Masks.MaskType == 'Dendrite')]['PathLength']))
        pixelSize = timeDict['pixelSize'][0]
        for center, dend in zip(DendCenters, DendPath):
            plt.text(center[0] / pixelSize + 18, center[1] / pixelSize, str(int(dend)), color='r')


def get_soma_dendrites(timeDict):
    if 'cellIndexAll' in timeDict:
        logger.info('Using cellIndexAll')
        soma_dends = np.unique(timeDict['Masks'].loc[timeDict['cellIndexAll']].DendNum.values)
    else:
        logger.info('Using cellIndex')
        soma_dends = np.unique(timeDict['Masks'].loc[timeDict['cellIndex']].DendNum.values)

    logger.info('Found %s as soma dendrite number(s)' % soma_dends)
    timeDict['soma_dends'] = soma_dends


def get_branch_id(timeDict):
    """

    :param timeDict:
    :return:
    """
    G = timeDict['graph']
    CellBodyPath = timeDict['CellBodyPath']
    xyStep = timeDict['UM_1X'] / timeDict['xyPixNum']
    Table2 = np.array(timeDict['InfoSWC']['TableSWC'])
    Table = Table2[:, 2:5]
    Table[:, 0] = Table[:, 0] * xyStep
    Table[:, 1] = Table[:, 1] * xyStep
    Table[:, 2] = Table[:, 2] * timeDict['anatomyZstep']
    GD = G.to_directed()
    dist = {}
    for node in GD:
        dist[node] = GD.number_of_edges(1, node)
    for node in GD:
        pred_nodes = [x for x in GD.predecessors(node)]
        for pred_node in pred_nodes:
            if CellBodyPath[pred_node] >= CellBodyPath[node]:
                GD.remove_edge(pred_node, node)
    dfs = nx.algorithms.traversal.dfs_preorder_nodes(GD, source=1.0)
    counterID = 1
    branchID = np.zeros(len(GD))
    for node in dfs:
        if node == 1:
            for nodeID in GD.successors(node):
                if nodeID != -1:
                    branchID[int(nodeID)] = counterID
                    counterID += 1
        else:
            succ = np.array(list(GD.successors(node)))
            if len(succ) == 1:
                branchID[int(succ[0])] = branchID[int(node)]
            elif len(succ) >= 2:
                out_deg = []
                for nodeID in succ:
                    GS = GD.subgraph(descendants(GD, nodeID))
                    ODS = GS.out_degree()
                    out_deg.append(np.sum([x[1] > 1 for x in list(ODS)]))
                max_deg = np.max(out_deg)
                putShaft = succ[(out_deg == max_deg).nonzero()[0]]
                if len(putShaft) == 1:
                    shaft_node = putShaft
                else:
                    angles = []
                    vec_in = Table[int(node) - 1, :] - Table[int(list(GD.predecessors(node))[0]) - 1, :]
                    for psn in putShaft:
                        vec_out = Table[int(psn) - 1, :] - Table[int(node) - 1, :]
                        angles.append(angle(vec_in, vec_out))
                    shaft_node = putShaft[np.argmin(angles)]
                for nodeID in succ:
                    if nodeID == shaft_node:
                        branchID[int(nodeID)] = branchID[int(node)]
                    else:
                        branchID[int(nodeID)] = counterID
                        counterID += 1
    branchID = branchID[1:]
    timeDict['branchID'] = branchID
    timeDict['directed_graph'] = GD


def get_branch_to_branch(timeDict, makeB=False):
    """ get connecting points, distance, relative branch order between branches

    :param timeDict:
    :return: dend_info
    dend_info axis 0: mask in dendrite 1, mask in dendrite 2, distance, soma traverse, relative branch order
    dend_info axis 1, 2: number of dendrites in session, the soma, number of terminal leafs
    """
    pathLenArray = timeDict['pathLengthAll']
    if makeB:
        excludeIndex = timeDict['excludeIndexB']
        Masks = timeDict['MasksB']
    else:
        excludeIndex = timeDict['excludeIndex']
        Masks = timeDict['Masks']
    if np.any(np.isnan(excludeIndex)):
        excludeIndex=[]
    G = timeDict['graph']
    branchID = timeDict['branchID']
    soma_dends = timeDict['soma_dends']
    exclude_dend = np.unique(Masks.loc[excludeIndex].DendNum.values)
    dend = timeDict['dend']
    num_dend = len(dend)
    transformInv = loadTransform(timeDict['path'])
    pixelSizeSession = timeDict['pixelSize']
    leafs = [x for x in G.nodes() if G.degree(x) == 1 and x != -1]
    n_total = num_dend + len(leafs)

    results = np.zeros((5, n_total, n_total))
    Table2 = np.array(timeDict['InfoSWC']['TableSWC'])
    Table = Table2[:, 2:5]
    Table[:, 0] = Table[:, 0]
    Table[:, 1] = Table[:, 1]
    Table[:, 2] = Table[:, 2]
    tree2 = KDTree(Table)
    swc_endpoints = Table2[np.array(leafs).astype(int) - 1, -1]
    index_list = []
    ids_list = []
    distance_list = []
    is_terminal = []
    branch_modal_node = []
    for i in range(num_dend):
        dend1 = copy.deepcopy(dend[i])
        dend1_fov = []
        if len(dend1.shape) < 2:
            dend1 = [dend1]
        for center in dend1:
            point = transformInv.TransformPoint(
                (float(center[2]), center[1] * pixelSizeSession[0], center[0] * pixelSizeSession[0]))
            dend1_fov.append((point[1], point[2], point[0]))
        distances, Indexes = tree2.query(dend1_fov)
        if np.median(distances) > 10:
            logger.info('Median dist 1 high: = %f' % np.median(distances))
        distances_dend1, Indexes_dend1 = tree2.query(dend1_fov)
        index_list.append(Indexes_dend1)
        distance_list.append(distances_dend1)
        Ids_dend1 = Table2[Indexes_dend1, 0]
        Ids_dend1 = get_excluded_index(Ids_dend1, pathLenArray, i)
        branch_ids, branch_ids_count = np.unique(branchID[np.array(Ids_dend1).astype(int) - 1], return_counts=True)
        high_id = branch_ids[np.argmax(branch_ids_count)]
        high_id2 = np.where(branchID[np.array(Ids_dend1).astype(int) - 1] == high_id)[0][0]
        branch_modal_node.append(Ids_dend1[high_id2])
        Ids_dend1 = np.unique(Ids_dend1)
        ids_list.append(Ids_dend1)
        is_terminal.append(len(np.intersect1d(Ids_dend1, swc_endpoints)) > 0)

    for x, y in itertools.product(range(n_total), range(n_total)):
        if x in soma_dends or y in soma_dends or x in exclude_dend or y in exclude_dend or x == y:
            continue
        #     print(x, y)
        if x < num_dend:
            distances_dend1 = distance_list[x]
            Indexes_dend1 = index_list[x]
            Ids_dend1 = ids_list[x]
            branch_modal_node_x = branch_modal_node[x]
        elif x == num_dend:
            Ids_dend1 = np.array([1])
            distances_dend1 = np.array([0])
            Indexes_dend1 = np.array([0])
            branch_modal_node_x = 1.0
        else:
            Ids_dend1 = np.array([leafs[x - num_dend]])
            distances_dend1 = np.array([0])
            Indexes_dend1 = np.array([leafs[x - num_dend] - 1]).astype(int)
            branch_modal_node_x = Ids_dend1[0]

        if y < num_dend:
            distances_dend2 = distance_list[y]
            Indexes_dend2 = index_list[y]
            Ids_dend2 = ids_list[y]
            branch_modal_node_y = branch_modal_node[y]
        elif y == num_dend:
            Ids_dend2 = np.array([1])
            distances_dend2 = np.array([0])
            Indexes_dend2 = np.array([0])
            branch_modal_node_y = 1.0
        else:
            Ids_dend2 = np.array([leafs[y - num_dend]])
            distances_dend2 = np.array([0])
            Indexes_dend2 = np.array([leafs[y - num_dend] - 1]).astype(int)
            branch_modal_node_y = Ids_dend2[0]

        dist = np.inf
        dend1_close = 0
        dend2_close = 0
        for x2, y2 in itertools.product(Ids_dend1, Ids_dend2):
            if timeDict['pathLengthAll'][int(x2 - 1), int(y2 - 1)] < dist:
                dist = timeDict['pathLengthAll'][int(x2 - 1), int(y2 - 1)]
                dend1_close = x2
                dend2_close = y2
        if np.isfinite(dist):
            Ids_dend1_best = np.where(Table2[Indexes_dend1, 0] == dend1_close)[0]
            dend1_best_dist = distances_dend1[Ids_dend1_best]
            dend1_best_offset = np.argmin(dend1_best_dist)
            dend1_best_point = Ids_dend1_best[dend1_best_offset]

            Ids_dend2_best = np.where(Table2[Indexes_dend2, 0] == dend2_close)[0]
            dend2_best_dist = distances_dend2[Ids_dend2_best]
            dend2_best_offset = np.argmin(dend2_best_dist)
            dend2_best_point = Ids_dend2_best[dend2_best_offset]

            node_list = nx.shortest_path(G, dend1_close, dend2_close)
            crossing_soma = 1.0 in node_list
            node_list2 = nx.shortest_path(G, branch_modal_node_x, branch_modal_node_y)
            relative_branch_order = len(np.unique(branchID[np.array(node_list2).astype(int) - 1])) - 1
            results[:, x, y] = np.array(
                [dend1_best_point, dend2_best_point, dist, crossing_soma, relative_branch_order])
        else:
            raise ValueError('Distance between %d and %d not finite' % (x, y))
    if makeB:
        timeDict['dend_infoB'] = results
        timeDict['dend_is_terminalB'] = is_terminal
    else:
        timeDict['dend_info'] = results
        timeDict['dend_is_terminal'] = is_terminal


def get_mask_to_mask(timeDict, makeB=False, smoothWin=50):
    """

    :param timeDict:
    :return:
    """

    All = []
    for x in timeDict['dend']:
        if len(x.shape) > 1:
            All.append(x[:, [1, 0, 2]])
        else:
            All.append(x[[1, 0, 2]])
    if makeB:
        Masks = timeDict['MasksB']
        excludeIndex = timeDict['excludeIndexB']
        results = timeDict['dend_infoB']
    else:
        Masks = timeDict['Masks']
        excludeIndex = timeDict['excludeIndex']
        results = timeDict['dend_info']

    G = timeDict['graph']
    AnatomyCenters = np.array([np.array(x) for x in Masks.Centers])
    # find closest point in interpolated space
    segList = []
    Dist = []
    Indexes = []
    maskOrder = []

    for segNum, segment in enumerate(All):
        maskIdx = (Masks.DendNum == segNum).values.nonzero()[0].astype(int)
        AnatomyCentersSeg = AnatomyCenters[maskIdx]
        if len(segment.shape) < 2:
            segment = segment[np.newaxis, :]
        x = segment[:, 0]
        y = segment[:, 1]
        z = segment[:, 2]
        # subtract one for matlab 1 indexed problem
        x = x - 1
        y = y - 1
        z = z - 1
        # move to real space coordinates (um)
        x = x * timeDict['pixelSize'][0]
        y = y * timeDict['pixelSize'][1]
        z = z * timeDict['pixelSize'][2]
        sx = meanSmooth(x, smoothWin)
        sy = meanSmooth(y, smoothWin)
        sz = meanSmooth(z, smoothWin)
        dist = np.hstack([0, np.cumsum(np.sqrt(np.diff(sx) ** 2 + np.diff(sy) ** 2 + np.diff(sz) ** 2))])
        tree = KDTree(np.array(list(zip(x, y, z))))
        distances, indexes = tree.query(AnatomyCentersSeg)
        Indexes = np.hstack([Indexes, len(Dist) + np.array(indexes)])
        Dist = np.hstack([Dist, dist])
        maskOrder = np.hstack([maskOrder, maskIdx])
        segList = np.hstack([segList, segNum * np.ones(len(dist))])

    result2 = copy.copy(results)
    Indexes = Indexes.astype(int)
    segList = segList.astype(int)
    maskOrder = maskOrder.astype(int)
    nDend = len(All)
    cumIdx = 0
    for i in range(result2.shape[1]):
        result2[0, i, :] = result2[0, i, :] + cumIdx
        result2[1, :, i] = result2[1, :, i] + cumIdx
        if i < nDend:
            if len(All[i].shape)>1:
                cumIdx = cumIdx + len(All[i])
            else:
                cumIdx = cumIdx + 1
    nMasks = len(Indexes)
    pathLengthGrid = np.ones([nMasks, nMasks]) * np.NAN
    somaTravGrid = np.ones([nMasks, nMasks]) * np.NAN
    relBranchOrder = np.ones([nMasks, nMasks]) * np.NAN
    for i in range(nMasks):
        for j in range(nMasks):
            segx = segList[Indexes[i]]
            segy = segList[Indexes[j]]
            if segx == segy:
                pathLengthGrid[maskOrder[i], maskOrder[j]] = np.absolute(Dist[Indexes[i]] - Dist[Indexes[j]])
            else:
                distx = np.absolute(Dist[Indexes[i]] - Dist[int(result2[0, segx, segy])])
                disty = np.absolute(Dist[Indexes[j]] - Dist[int(result2[1, segx, segy])])
                pathLengthGrid[maskOrder[i], maskOrder[j]] = distx + disty + result2[2, segx, segy]
            somaTravGrid[maskOrder[i], maskOrder[j]] = result2[3, segx, segy]
            relBranchOrder[maskOrder[i], maskOrder[j]] = result2[4, segx, segy]
    pathLengthSoma = np.ones(nMasks) * np.NAN
    branchOrder = np.ones(nMasks) * np.NAN
    for i in range(nMasks):
        segx = segList[Indexes[i]]
        distx = np.absolute(Dist[Indexes[i]] - Dist[int(result2[0, segx, nDend])])
        pathLengthSoma[maskOrder[i]] = distx + result2[2, segx, nDend]
        branchOrder[maskOrder[i]] = result2[4, segx, nDend] - 1

    leafs = [x for x in G.nodes() if G.degree(x) == 1 and x != -1]
    nTerminal = len(leafs)
    pathLengthTerminals = np.ones([nMasks, nTerminal]) * np.NAN
    pathLengthTerminalsTrav = np.ones([nMasks, nTerminal]) * np.NAN
    for i in range(nMasks):
        for j in range(nTerminal):
            segx = segList[Indexes[i]]
            distx = np.absolute(Dist[Indexes[i]] - Dist[int(result2[0, segx, nDend + j])])
            pathLengthTerminals[maskOrder[i], j] = distx + result2[2, segx, nDend + j]
            pathLengthTerminalsTrav[maskOrder[i], j] = result2[3, segx, nDend + j]
    if np.any(np.isfinite(excludeIndex)):
        pathLengthGrid[:, excludeIndex] = np.NAN
        pathLengthGrid[excludeIndex, :] = np.NAN
        pathLengthSoma[excludeIndex] = np.NAN
        branchOrder[excludeIndex] = np.NAN
        pathLengthTerminals[excludeIndex, :] = np.NAN
        pathLengthTerminalsTrav[excludeIndex, :] = np.NAN
        somaTravGrid[:, excludeIndex] = np.NAN
        somaTravGrid[excludeIndex, :] = np.NAN
        relBranchOrder[:, excludeIndex] = np.NAN
        relBranchOrder[excludeIndex, :] = np.NAN

    if makeB:
        timeDict['pathLengthGridB'] = pathLengthGrid
        timeDict['pathLengthSomaB'] = pathLengthSoma
        timeDict['branchOrderB'] = branchOrder
        timeDict['pathLengthTerminalsB'] = pathLengthTerminals
        timeDict['pathLengthTerminalsTravB'] = pathLengthTerminalsTrav
        timeDict['somaTravGridB'] = somaTravGrid
        timeDict['relBranchOrderB'] = relBranchOrder
    else:
        timeDict['pathLengthGrid'] = pathLengthGrid
        timeDict['pathLengthSoma'] = pathLengthSoma
        timeDict['branchOrder'] = branchOrder
        timeDict['pathLengthTerminals'] = pathLengthTerminals
        timeDict['pathLengthTerminalsTrav'] = pathLengthTerminalsTrav
        timeDict['somaTravGrid'] = somaTravGrid
        timeDict['relBranchOrder'] = relBranchOrder


def findPathLength(timeDict, session, makeB=False, plot=True, smoothWin=50):
    """ calculates the path length along the dendrite to all masks 0 is cell body

    :param timeDict: time course dictionary
    :param session:
    :param loadSWC: load data from swc reconstruction file
    :return:
    """
    get_graph_from_swc(timeDict, session)
    get_path_length(timeDict, makeB=makeB, plot=plot)
    get_branch_id(timeDict)
    get_soma_dendrites(timeDict)
    get_branch_to_branch(timeDict, makeB=makeB)
    get_mask_to_mask(timeDict, makeB=makeB, smoothWin=smoothWin)


def meanSmooth(x, winSize):
    halfWin = int(winSize / 2)
    X = np.ones([len(x), winSize + 1]) * np.NAN
    for i in range(halfWin):
        X[(halfWin - i):, i] = x[:-(halfWin - i)]
        X[-(halfWin - i):, i] = np.NAN
    X[:, halfWin] = x
    for i in range(1, halfWin + 1):
        X[:-i, i + halfWin] = x[i:]
        X[:i, i + halfWin] = np.NAN
    return np.nanmean(X, axis=1)


def get_excluded_index(dendrite_ids, path_length, dendrite_number, dist_threshold=10):
    """ checks for discontinuities in distance between reconstruction points along a dendrites and tries to solve them
    by either cutting from the first one onwards or cutting from the beginning

    :param dendrite_ids: reconstruction points for current dendrite
    :param path_length: path length along the dendrite between reconstruction points
    :param dendrite_number: dendrite number in session
    :param dist_threshold: two reconstruction points with a distance larger then dist count as a discontinuity
    :return: reconstruction points after deleting the discontinuities
    """
    # get distances between points
    dist_dend = []
    for i in range(len(dendrite_ids) - 1):
        dist_dend.append(path_length[int(dendrite_ids[i]-1), int(dendrite_ids[i + 1]-1)])
    large_dist = (np.array(dist_dend) > dist_threshold).nonzero()[0]
    n_pairs = len(large_dist) / 2
    # try cutting from the first one
    exclude_index_1 = []
    for i in range(int(n_pairs)):
        exclude_index_1.extend(range(large_dist[i * 2] + 1, large_dist[i * 2 + 1] + 1))
    if len(large_dist) % 2:
        exclude_index_1.extend(range(large_dist[-1] + 1, len(dendrite_ids) + 1))
    large_dist2 = np.hstack((0, large_dist))
    n_pairs2 = len(large_dist2) / 2
    # try cutting from the beginning
    exclude_index_2 = []
    for i in range(int(n_pairs2)):
        exclude_index_2.extend(range(large_dist2[i * 2] + 1, large_dist2[i * 2 + 1] + 1))
    if len(large_dist2) % 2:
        exclude_index_2.extend(range(large_dist2[-1] + 1, len(dendrite_ids) + 1))

    # The better one is the longer one (smaller exclusion points)
    if len(exclude_index_1) > 0 or len(exclude_index_2) > 0:
        if len(exclude_index_1) < len(exclude_index_2):
            if len(exclude_index_1) > 0:
                logger.info('Excluded dend Dist 1 %d: %s' % (dendrite_number, exclude_index_1))
                dendrite_ids = np.delete(dendrite_ids, exclude_index_1)
        elif len(exclude_index_2) > 0:
            logger.info('Excluded dend Dist 2 %d: %s' % (dendrite_number, exclude_index_2))
            dendrite_ids = np.delete(dendrite_ids, exclude_index_2)
    return dendrite_ids


def showDendrites(timeDict, showB=False, exclude=None):
    """ helper function to plot the dendrite numbers

    :param timeDict: time course dictionary
    """
    if not showB:
        Masks = timeDict['Masks']
        labelimgAll = timeDict['labelimgAll']
        dendLabelTable = timeDict['dendLabelTable']
    else:
        Masks = timeDict['MasksB']
        labelimgAll = timeDict['labelimgAllB']
        dendLabelTable = timeDict['dendLabelTableB']
    plt.figure(figsize=(12, 12))
    plt.imshow(labelimgAll.max(axis=2).transpose(1, 0))
    DendCenters = np.asarray(list(Masks[(Masks.MaskType == 'Dendrite')]['Centers']))
    DendPath = np.asarray(list(Masks[(Masks.MaskType == 'Dendrite')].index))
    pixelSize = timeDict['pixelSize'][0]
    for center, dend in zip(DendCenters, DendPath):
        plt.text(center[0] / pixelSize + 18, center[1] / pixelSize, str(int(dend)), color='r')
    dendNums = np.unique(dendLabelTable[:, 1].T)
    for num in dendNums:
        if exclude is None or num not in exclude:
            current = Masks[(Masks.MaskType == 'Dendrite') & (Masks.DendNum == num)].Centers
            XY = (np.array(list(map(list, current.values)))[:, [0, 1]].mean(axis=0) / pixelSize).astype(int)
            plt.text(XY[0], XY[1], str(num), size=20, color='white')


def showSpines(timeDict, selected=None, pixSize=376.666, showB=False):
    """

    :param timeDict:
    :return:
    """
    if showB:
        masks = timeDict['MasksB']
    else:
        masks = timeDict['Masks']
    shapes = []
    annotations = []
    for mark in selected:
        x = masks.loc[mark].Centers[0] * 1000.0 / pixSize
        y = masks.loc[mark].Centers[1] * 1000.0 / pixSize
        shapes.append({
            'type': 'circle',
            'xref': 'x',
            'yref': 'y',
            'x0': x - 4,
            'y0': y - 4,
            'x1': x + 4,
            'y1': y + 4,
            'line': {
                'color': 'rgba(50, 171, 96, 1)',
            },
        })
        annotations.append(dict(
            x=x,
            y=y + 5,
            xref='x',
            yref='y',
            text=str(mark),
            font=dict(
                family='Courier New, monospace',
                size=16,
                color='#ffffff'),
            showarrow=True,
            arrowhead=7,
            ax=0,
            ay=-40
        ))
    if showB:
        labelimgAll = timeDict['labelimgAllB']
    else:
        labelimgAll = timeDict['labelimgAll']
    layout = go.Layout(height=labelimgAll.shape[1] * 2,
                       width=labelimgAll.shape[0] * 2,
                       shapes=shapes,
                       annotations=annotations)
    heatmap = labelimgAll.max(axis=2).transpose(1, 0)
    trace = go.Heatmap(z=heatmap)
    py.iplot(dict(data=[trace], layout=layout))


def get_path_length_grid(timeDict):
    """

    :param timeDict:
    :return:
    """
    Masks = timeDict['Masks']
    if isinstance(timeDict['pathLengthAll'], dict):
        a = np.zeros((len(timeDict['pathLengthAll']) + 1, len(timeDict['pathLengthAll']) + 1))
        for key, value in timeDict['pathLengthAll'].iteritems():
            for key2 in sorted(value):
                a[int(key), int(key2)] = value[key2]
        timeDict['pathLengthAll'] = a
    pathLength = timeDict['pathLengthAll']
    pathLengthGrid = np.zeros((len(Masks), len(Masks)))
    for x in range(len(Masks)):
        xNode = int(Masks.loc[x].ParentId)
        for y in range(len(Masks)):
            yNode = int(Masks.loc[y].ParentId)
            if x != y:
                pathLengthGrid[x, y] = pathLength[xNode - 1, yNode - 1]
    excludeIndex = timeDict['excludeIndex']
    somaPathLength = np.array(timeDict['Masks']['PathLength'])
    spineIdx = np.array(timeDict['Masks']['MaskType'] == 'Spine').nonzero()[0]
    somaPathLength[spineIdx] = somaPathLength[spineIdx]
    if np.any(np.isfinite(excludeIndex)):
        pathLengthGrid[:, excludeIndex] = np.NAN
        # pathLengthGrid[excludeIndex, :] = np.NAN
    timeDict['pathLengthGrid'] = pathLengthGrid


def get_high_res_path_length_grid(timeDict, smoothWin=50):
    def meanSmooth(x,winSize):
        X = np.ones([len(x),winSize+1])*np.NAN
        for i in range(winSize/2):
            X[((winSize/2)-i):,i]=x[:-((winSize/2)-i)]
            X[-(winSize/2-i):,i]=np.NAN
        X[:,winSize/2]=x
        for i in range(1,winSize/2+1):
            X[:-i,i+winSize/2]=x[i:]
            X[:i,i+winSize/2]=np.NAN
        return np.nanmean(X,axis=1)

    All=[]
    for x in timeDict['dend']:
        if len(x.shape)>1:
            All.append(x[:,[1,0,2]])
        else:
            All.append(x[[1,0,2]])
    Masks = timeDict['Masks']
    nMasks = len(Masks)
    AnatomyCenters = np.array([np.array(x) for x in Masks.Centers])
    # find closest point in interpolated space
    segList = np.zeros(nMasks)
    Dist = np.zeros(nMasks)

    for segNum, segment in enumerate(All):
        maskIdx = (Masks.DendNum==segNum).values.nonzero()[0]
        AnatomyCentersSeg = AnatomyCenters[maskIdx]
        if len(segment.shape)<2:
            Dist[maskIdx]=np.zeros(len(maskIdx))
            segList[maskIdx]=segNum * np.ones(len(maskIdx))
            continue
        x = segment[:,0]
        y = segment[:,1]
        z = segment[:,2]
        # subtract one for matlab 1 indexed problem
        x = x - 1
        y = y - 1
        z = z - 1
        # move to real space coordinates (um)
        x = x * timeDict['pixelSize'][0]
        y = y * timeDict['pixelSize'][1]
        z = z * timeDict['pixelSize'][2]
        sx = meanSmooth(x,smoothWin)
        sy = meanSmooth(y,smoothWin)
        sz = meanSmooth(z,smoothWin)
        dist = np.hstack([0,np.cumsum(np.sqrt(np.diff(sx)**2+np.diff(sy)**2+np.diff(sz)**2))])
        tree = KDTree(np.array(list(zip(x, y, z))))
        distances, indexes = tree.query(AnatomyCentersSeg)
        Dist[maskIdx]=dist[indexes]
        segList[maskIdx]=segNum * np.ones(len(maskIdx))
    pathLengthGrid=np.ones([nMasks,nMasks])*np.NAN
    for i in range(nMasks):
        for j in range(nMasks):
            if segList[i]==segList[j]:
                pathLengthGrid[i,j]=np.absolute(Dist[i]-Dist[j])
    excludeIndex = timeDict['excludeIndex']
    if np.any(np.isfinite(excludeIndex)):
        pathLengthGrid[:, excludeIndex] = np.NAN
    timeDict['pathLengthGridHR'] = pathLengthGrid


def get_branch_dist(timeDict, plot=False, warn_dist=10):
    # Note: Not working most of the time
    dend = timeDict['dend']
    num_dend = len(dend)
    transformInv = loadTransform(timeDict['path'])
    pixelSizeSession = timeDict['pixelSize']
    xyStep = timeDict['UM_1X'] / timeDict['xyPixNum']

    # find closest point in interpolated space
    img = tf.imread(timeDict['SessionPath']).astype(int).transpose(1, 2, 0)
    imgSize = img.shape
    Info = timeDict['InfoSWC']
    All = Info['AllSWC']
    Connecting = list([])
    segList = list([])
    for segNum, segment in enumerate(All):
        x, y, z = np.unravel_index(segment, imgSize, 'F')
        # subtract one for matlab 1 indexed problem
        x = x - 1
        y = y - 1
        z = z - 1
        # move to real space coordinates (um)
        x = x * xyStep
        y = y * xyStep
        z = z * timeDict['anatomyZstep']
        if len(x.shape) > 0:
            Connecting.extend(list(zip(x, y, z)))
            segList.extend([segNum] * x.shape[0])
        else:
            Connecting.append((x, y, z))
            segList.append(segNum)
    Connecting = np.array(Connecting)
    tree = KDTree(Connecting)

    Table2 = np.array(timeDict['InfoSWC']['TableSWC'])
    Table = Table2[:, 2:5]
    Table[:, 0] = Table[:, 0]  # * xyStep
    Table[:, 1] = Table[:, 1]  # * xyStep
    Table[:, 2] = Table[:, 2]  # * timeDict['anatomyZstep']
    tree2 = KDTree(Table)
    # results in ((point in dend x, point in dend y, dist), dend x, dend y)
    results = np.zeros((3, num_dend, num_dend))
    for x, y in itertools.product(range(num_dend), range(num_dend)):
        dend1 = copy.deepcopy(dend[x]).astype('float')
        dend2 = copy.deepcopy(dend[y]).astype('float')
        dend1_fov = []
        dend2_fov = []
        for center in dend1:
            point = transformInv.TransformPoint(
                (center[2], center[1] * pixelSizeSession[1], center[0] * pixelSizeSession[0]))
            dend1_fov.append((point[1], point[2], point[0]))
        for center in dend2:
            point = transformInv.TransformPoint(
                (center[2], center[1] * pixelSizeSession[1], center[0] * pixelSizeSession[0]))
            dend2_fov.append((point[1], point[2], point[0]))
        distances, indexes = tree.query(dend1_fov)
        if np.median(distances) > warn_dist:
            logger.info('Median dist high: = %f' % np.median(distances))
        distances_dend1, indexes_dend1 = tree2.query(dend1_fov)
        Ids_dend1 = Table2[indexes_dend1, 0]
        distances_dend2, indexes_dend2 = tree2.query(dend2_fov)
        Ids_dend2 = Table2[indexes_dend2, 0]
        dendNum_dend2 = Table2[indexes_dend2, -1]
        dendNum_dend1 = Table2[indexes_dend1, -1]
        test = copy.deepcopy(dendNum_dend1)
        uD, counts = np.unique(test, return_counts=True)
        exIdx = []
        for k in uD:
            check = (test == k).nonzero()[0]
            if len(check) > 1:
                if np.max(np.diff(check)) > 1:
                    exIdx = (test != uD[np.argmax(counts)]).nonzero()[0]
            else:
                exIdx = (test != uD[np.argmax(counts)]).nonzero()[0]
        if len(exIdx) > 0:
            logger.info('Excluded dend1: %s' % exIdx)
            Ids_dend1 = np.delete(Ids_dend1, exIdx)
        test = copy.deepcopy(dendNum_dend2)
        uD, counts = np.unique(test, return_counts=True)
        exIdx = []
        for k in uD:
            check = (test == k).nonzero()[0]
            if len(check) > 1:
                if np.max(np.diff(check)) > 1:
                    exIdx = (test != uD[np.argmax(counts)]).nonzero()[0]
            else:
                exIdx = (test != uD[np.argmax(counts)]).nonzero()[0]
        if len(exIdx) > 0:
            logger.info('Excluded dend2: %s' % exIdx)
            Ids_dend2 = np.delete(Ids_dend2, exIdx)

        Ids_dend2 = np.unique(Ids_dend2).astype('int')
        Ids_dend1 = np.unique(Ids_dend1).astype('int')
        dist = np.inf
        dend1_close = 0
        dend2_close = 0
        for x2, y2 in itertools.product(Ids_dend1, Ids_dend2):
            if timeDict['pathLengthAll'][x2][y2] < dist:
                dist = timeDict['pathLengthAll'][x2][y2]
                dend1_close = x2
                dend2_close = y2
        if np.isfinite(dist):
            Ids_dend1_best = np.where(Table2[indexes_dend1, 0] == dend1_close)[0]
            dend1_best_dist = distances_dend1[Ids_dend1_best]
            dend1_best_offset = np.argmin(dend1_best_dist)
            dend1_best_point = Ids_dend1_best[dend1_best_offset]

            Ids_dend2_best = np.where(Table2[indexes_dend2, 0] == dend2_close)[0]
            dend2_best_dist = distances_dend2[Ids_dend2_best]
            dend2_best_offset = np.argmin(dend2_best_dist)
            dend2_best_point = Ids_dend2_best[dend2_best_offset]
            results[:, x, y] = np.array([dend1_best_point, dend2_best_point, dist])
        else:
            results[:, x, y] = np.array([np.nan, np.nan, dist])
    if plot:
        plt.imshow(results[2, :, :])
        plt.colorbar(label='Distance')
        plt.xlabel('Dend#')
        plt.ylabel('Dend#')
    return results


def get_mask_index(timeDict, mask='Spine', use_B=False, noise_th=None):
    """

    :param timeDict: timeDict to use
    :param mask: options are 'Spine' and 'Dendrite'
    :param use_B: Make masksB etc.
    :param noise_th: if None will return all mask index if float will return mean noise < then threshold
    :return: index of masks
    """
    if use_B:
        b = 'B'
    else:
        b = ''
    masks = timeDict['Masks' + b]
    exclude = timeDict['excludeIndex' + b]
    indexs_all = np.where(masks.MaskType == mask)[0]
    indexs_good = np.setdiff1d(indexs_all, exclude)
    if noise_th is not None:
        noise = np.nanmean(timeDict['TCNoise' + b], axis=1)
        good_noise = np.where(noise < noise_th)[0]
        return np.intersect1d(indexs_good, good_noise)
    else:
        return indexs_good
