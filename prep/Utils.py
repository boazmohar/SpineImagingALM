# Encoding: utf-8
""" General utilities for Thunder Images objects and misc.

"""
from __future__ import print_function

import copy
import io
import glob
import logging
import math
import os
import re
import sys
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import requests
import pandas as pd
import pyspark
import thunder as td
import matplotlib.patches as mpatches
from IPython.display import display as ipython_display
from future.utils import iteritems
from ipywidgets import widgets
from ipywidgets.widgets import interactive, interact
from matplotlib.offsetbox import AnchoredOffsetbox
from pylab import rms_flat
from skimage.external import tifffile
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from skimage.measure import label
from skimage.restoration.inpaint import inpaint_biharmonic
from pySparkUtils.utils import thunder_decorator


try:
    basestring = basestring
except NameError:
    basestring = str

# Setup logging
logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)

try:
    from ScanImageTiffReader import ScanImageTiffReader
except OSError:
    logger.error('wrong ScanImageTiffReader DLL?')
    ScanImageTiffReader = None


@thunder_decorator
def records_per_partitions(rdd):
    """

    :param obj: Spark rdd
    :return: number of records per partition
    """
    return rdd.glom().map(len).collect()


def anova_from_stats(m, sd, n):
    ##(a) anova table
    k = float(len(m))  # number of groups
    Xg = np.sum(n * m) / np.sum(n)

    dfb = k - 1  # degree of freedom
    dfw = np.sum(n) - k  # degree of freedom

    MSb = np.sum(n * (m - Xg) ** 2) / (k - 1)  # MS between
    MSw = np.sum((n - 1) * sd ** 2) / dfw  # MS within
    SSb = dfb * MSb
    SSw = dfw * MSw
    SSt = SSb + SSw

    return MSb / MSw  # f value


def get(object_):
    what(object_)
    return object_


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def get_joint_valid(x1, y1):
    x_valid = np.logical_not(np.isnan(x1)).nonzero()[0]
    y_valid = np.logical_not(np.isnan(y1)).nonzero()[0]
    return np.intersect1d(x_valid, y_valid)


def select_session(base='/groups/svoboda/svobodalab/users/Aaron', session_class='Spine'):
    """ plots a widget that enables selection of session from a base directory
    The convention is:
    1. folder with date_animal (170506_BMWR56)
    2. folder with Run number inside contains the raw data(170506_BMWR56/Run1, 160214_BMWR26/Stack2)
    3. pickle file with in the animal folder which starts with the run (170506_BMWR56/Run1Step3.p)
    
    :param base: where to look for sessions
    :param session_class: which type of session object to load
    :return: the button (load.session will be the session object after pressing load)
    """
    if session_class == 'Spine':
        from prep.SpineSession import SpineSession as Session
    elif session_class == 'Base':
        from prep.Session import Session as Session
    elif session_class == 'Vis':
        from prep.VisSession import VisSession as Session

    def list_dirs(path):
        return [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]

    def list_sessions(animal):
        dates.options = sorted([x[1] for x in dir_list if x[2] == animal])
        if dates.options is not None:
            dates.value = dates.options[0]
        list_runs(dates.options[0])

    def list_runs(date):
        animal = animals.value
        if animal is not None:
            dir1 = [x[0] for x in dir_list if x[1] == date and x[2] == animal]
            if dir1 is not None:
                dir1 = dir1[0]
            runs.options = sorted(list_dirs(os.path.join(base, dir1)))

    def list_files(run):
        animal = animals.value
        date = dates.value
        dir1 = [x[0] for x in dir_list if x[1] == date and x[2] == animal]
        if dir1 is not None:
            dir1 = dir1[0]
            s = os.path.join(base, dir1, '')
            files.options = sorted(map(lambda y: os.path.basename(y)[len(run):-2], glob.glob(s + run + '*.p')))
        files.options = ['None'] + list(files.options)

    dirs = list_dirs(base)
    r = re.compile(r"([0-9]{6})_(.*)")
    dir_list = [(m.group(0), m.group(1), m.group(2)) for l in dirs for m in [r.search(l)] if m]
    animal_list = sorted(list(set(map(lambda x: x[2], dir_list))))
    animals = widgets.Dropdown(options=animal_list)
    wA = interactive(list_sessions, animal=animals, button_style='danger')
    dates = widgets.Dropdown(options=[])
    wS = interactive(list_runs, date=dates)
    runs = widgets.Dropdown(options=[], )
    wR = interactive(list_files, run=runs)
    files = widgets.Dropdown(options=[], )
    wF = interactive(lambda my_file: my_file, my_file=files)
    load = widgets.Button(description='Load', button_style='danger')
    ipython_display(wA, wS, wR, wF, load)

    def on_button_clicked(b):
        animal = wA.kwargs['animal']
        date = wS.kwargs['date']
        run = wR.kwargs['run']
        file_name = wF.kwargs['my_file']
        b.session = Session(basePath=base, animalID=animal, date=date, run=run)
        if file_name is not None and not file_name == 'None':
            b.session = b.session.load(file_name)

    load.on_click(on_button_clicked)
    return load


def convert_8bit(image, sat_percent=0.1, ignore_nans=True, ignore_zero=False):
    """

    :param image: input image to convert
    :param sat_percent: percent saturation
    :param ignore_nans: if True will calculate the percentage on non-nan pixels
    :return: image in uint8 and after saturation
    """
    if ignore_zero:
        image = copy.deepcopy(image)
        image[image == 0.0] = np.nan
    if ignore_nans:
        n_pix = np.sum(np.isfinite(image.flatten()))
        n_pix_sat = int(sat_percent / 100.0 * n_pix)
        sat_p = (1.0 - n_pix_sat / np.prod(image.shape)) * 100.0
    else:
        sat_p = 100.0 - sat_percent
    logger.info('Saturation set at %.4f percent' % sat_p)
    image_8bit = copy.deepcopy(image)
    image_8bit = image_8bit - np.nanmin(image_8bit)
    image_8bit = image_8bit / np.nanpercentile(image_8bit, sat_p)
    image_8bit[image_8bit > 1.0] = 1.0
    image_8bit = image_8bit * 255.0
    return image_8bit.astype(np.uint8)


def getLabelImg(sc, session, start, stop, TC, clipTC=(-2, 20), subMeansFlag=False, useC=True, t_index=None, useB=False):
    def labelTC(key):
        out = np.zeros(labelimgAllBC.value.shape)
        out[:] = np.NAN
        for j in range(len(cIdx)):
            out[labelimgAllBC.value == (cIdx[j] + 1)] = ZtBC.value[j, key]
        return key, np.nanmax(out, axis=2)

    def subMeans(Z):
        Ztsub = np.zeros(Z.shape)
        z_mu = np.nanmean(Z, axis=1)
        for j in range(Z.shape[0]):
            if not np.isnan(z_mu[j]):
                Ztsub[j, :] = Z[j, :] - z_mu[j]
        z_mu = np.nanmean(Z, axis=0)
        for j in range(Z.shape[1]):
            if not np.isnan(z_mu[j]):
                Ztsub[:, j] = Z[:, j] - z_mu[j]
        return Ztsub

    timeDict = session.timeDict
    if useC:
        c = copy.deepcopy(timeDict['c'])
    else:
        c = range(TC.shape[0])

    TC[TC < clipTC[0]] = clipTC[0]
    TC[TC > clipTC[1]] = clipTC[1]
    Zt0 = TC[c, :]
    if subMeansFlag:
        Zt0 = subMeans(Zt0)
    ZtBC = sc.broadcast(Zt0)
    if useB:
        labelimgAllBC = sc.broadcast(timeDict['labelimgAllB'])
    else:
        labelimgAllBC = sc.broadcast(timeDict['labelimgAll'])
    cIdx = c
    if t_index is None:
        dispRange = range(start, stop)
        tIdx = dispRange
    else:
        tIdx = t_index
    # idxList = zip(range(len(tIdx)))
    labelTCDict = sc.parallelize(tIdx, len(tIdx)).map(labelTC).collectAsMap()
    labelTCdiv = np.zeros((len(tIdx), labelimgAllBC.value.shape[0], labelimgAllBC.value.shape[1]))
    for j, i in enumerate(labelTCDict.keys()):
        labelTCdiv[j, :, :] = labelTCDict[i]
    return labelTCdiv


def adjust_spines(ax=None, spines=('bottom', 'left'), separation=10, smart_bounds=True):
    """
    
    :param ax: axis to adjust. If None will use plt.gca()
    :param spines: which spines out of: ['bottom','top','left','right']
    :param separation: move out y axis in points
    :param smart_bounds: apply smart bounds
    """
    if ax is None:
        ax = plt.gca()
    for loc, spine in ax.spines.items():
        if loc in spines:
            if separation > 0:
                spine.set_position(('outward', separation))  # outward by 10 points
            spine.set_smart_bounds(smart_bounds)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no y axis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no x axis ticks
        ax.xaxis.set_ticks([])


def what(obj):
    print('type: %s' % type(obj))
    if hasattr(obj, '__len__'):
        print('len: %s' % len(obj))
    if isinstance(obj, dict):
        print('key: value')
        print('~~~~~~~~~~')
        for key, value in iteritems(obj):
            print('%s: %s' % (key, type(value)), end=', ')
            if isinstance(value, np.ndarray):
                print(', shape: %s, size: %.2fMB' % (value.shape, value.nbytes * 1e-6))
            elif isinstance(value, (list, tuple)):
                print('len: %s: type: %s' % (len(value), type(value[0])))
            else:
                print('')
    else:
        if isinstance(obj, np.ndarray):
            print('shape: %s, dtype: %s, size: %.2fMB' % (str(obj.shape), str(obj.dtype), obj.nbytes * 1e-6))
        elif isinstance(obj, (str, basestring)):
            print(str)


def nanMeanByIndex(data, index, largeGroupMean=True):
    """ try to aggregate similar time points and nanmean then

    :param data: images object
    :param index: array of size (1, data.shape[0]) to aggregate by
    :param largeGroupMean: if True wil just mean. if False will run a local target type function (see getTarget)
    :return: numpy array of size (unique(index), data.shape[1], ...)
    """
    th = int((3.2e10 / 32) / np.prod(data.shape[1:])) // 2
    data = data.tordd()
    index = index.astype(int)

    def map_keys(x):
        if isinstance(x[0], tuple):
            x = (int(x[0][0]),  x[1])
        y = ((index[int(x[0])],), x[1])
        return y

    data = data.map(map_keys)

    def comb1(x, y):
        if len(x.shape) == 1:
            return y
        elif len(y.shape) == 1:
            return x
        elif len(x.shape) == 3 and len(y.shape) == 3:
            return np.stack([x, y])
        else:
            if len(x.shape) < len(y.shape):
                x = np.expand_dims(x, 0)
            elif len(y.shape) < len(x.shape):
                y = np.expand_dims(y, 0)
            z = np.concatenate((x, y))
            # if len(z.shape) == 4 and z.shape[0] > 100:
            #     if largeGroupMean:
            #         z = np.expand_dims(np.nanmean(z, axis=0), 0)
            return z

    def getTargetLocal(array):
        sz = array.shape
        if len(sz) == 3:
            return array
        result = np.zeros(shape=(sz[1], sz[2], sz[3]), dtype=array.dtype)
        Mean = np.nan_to_num(copy.deepcopy(
            np.nanmean(array, axis=0).reshape(sz[1] * sz[2], sz[3])))
        array2 = np.nan_to_num(copy.deepcopy(array.reshape(sz[0], sz[1] * sz[2], sz[3])))
        for i in range(sz[3]):
            CC = list()
            for k in range(sz[0]):
                CC.append(np.corrcoef(array2[k, :, i], Mean[:, i])[0, 1])
            CC = np.array(CC)
            if sz[0] < 30:
                points = sz[0]
            elif sz[0] < 100:
                points = np.round(sz[0] / 2).astype(int)
            elif sz[0] < 200:
                points = np.round(sz[0] / 3).astype(int)
            elif sz[0] < 300:
                points = np.round(sz[0] / 4).astype(int)
            else:
                points = np.round(sz[0] / 5).astype(int)
            ind = np.argpartition(CC, -points)[-points:]
            result[:, :, i] = np.nanmean(array[ind, :, :, i], axis=0).astype('float32')
        return result

    def getMean(array):
        sz = array.shape
        if len(sz) == 3:
            return array
        else:
            return np.nanmean(array, axis=0)

    # If the data (i.e. a single group) is bigger then 4GB (3.2e10 bits) the aggregation will fail in spark,
    # so split it into two or more group and average the result

    index2 = copy.deepcopy(index)
    counts = np.bincount(index2)
    bigGroups = np.where(counts > th)[0]
    fixList = []
    for bigGroup in bigGroups:
        index2 = np.where(index == bigGroup)[0]
        numGroups = len(index2) // th
        for k in range(1, numGroups):
            newVal = np.max(index) + 1
            index[index2[k::numGroups]] = newVal
            fixList.append((bigGroup, copy.deepcopy(index2[k::numGroups])))

    data = data.aggregateByKey(np.array([]), comb1, comb1)
    if largeGroupMean:
        data = data.mapValues(getMean).collectAsMap()
    else:
        data = data.mapValues(getTargetLocal).collectAsMap()
    r = np.array([data[idx] for idx in sorted(data.keys())])
    extraIndex = r.shape[0] - len(fixList)
    extra = r[extraIndex:, ...]
    r = r[:extraIndex, ...]
    for k, (bigGroup, index2) in enumerate(fixList):
        comb = np.nanmean(np.stack((r[bigGroup, ...], extra[k, ...]), 0), 0)
        r[bigGroup, ...] = comb
    return r


def getCrop(data, index=None, yStart=2, yStop=None, xStart=0, xStop=None, tIndex=None, display=0.5, excludeIndex=True):
    """ shows the data after a crop to x, y, z

    Parameters
    ----------
    data: volume to use
    index: index of z planes to exclude
    yStart: axis 1 start
    yStop: axis 1 stop
    xStart: axis 2 start
    xStop: axis 3 stop
    tIndex: indexes in time to exclude
    display: saturation of display, between (0, 1] lower is more saturated
    excludeIndex: if true index planes are excluded from the crop

    Returns
    -------
    A dictionary of the resulting crop
    """
    sz = data.shape
    if yStop is None:
        yStop = sz[0]
    if xStop is None:
        xStop = sz[1]
    ySlice = slice(yStart, yStop)
    xSlice = slice(xStart, xStop)
    zIndex = np.ones(sz[2], dtype=bool)
    if index is not None:
        zIndex[index] = False
    if not excludeIndex:
        zIndex = np.logical_not(zIndex)
    if display is not None:
        data2 = copy.deepcopy(data[ySlice, xSlice, zIndex])
        numPlots = np.ceil(data2.shape[2] ** 0.5).astype(int)
        Min = np.nanmin(data2)
        Max = np.nanmax(data2) * display
        if numPlots > 1:
            fig, axes = plt.subplots(numPlots, numPlots, figsize=(16, 8), subplot_kw={'xticks': [], 'yticks': []})
            fig.subplots_adjust(hspace=0.05, wspace=0.05)
            for i, ax in enumerate(axes.flat):
                if i < data2.shape[2]:
                    ax.imshow(data2[:, :, i], vmin=Min, vmax=Max, cmap='gray', interpolation='nearest')
                    ax.set_title(np.where(zIndex)[0][i], color='white', y=0.8, x=0.05)
        else:
            plt.imshow(data2[:, :, 0], vmin=Min, vmax=Max, cmap='gray', interpolation='nearest')

    crop = dict()
    crop['tIndex'] = tIndex
    crop['zIndex'] = zIndex
    crop['ySlice'] = ySlice
    crop['xSlice'] = xSlice
    return crop


def xCorrelation(x, y=None, maxlags=None, norm='biased'):
    """Cross-correlation using numpy.correlate

    Estimates the cross-correlation (and autocorrelation) sequence of a random
    process of length N. By default, there is no normalisation and the output
    sequence of the cross-correlation has a length 2*N+1.

    :param array x: first data array of length N
    :param array y: second data array of length N. If not specified, computes the
        autocorrelation.
    :param maxlags: compute cross correlation between [-maxlags:maxlags]
        when maxlags is not specified, the range of lags is [-N+1:N-1].
    :param norm: normalisation in ['biased', 'unbiased', None, 'coeff']

    The true cross-correlation sequence is

    .. math:: r_{xy}[m] = E(x[n+m].y^*[n]) = E(x[n].y^*[n-m])

    However, in practice, only a finite segment of one realization of the
    infinite-length random process is available.

    The correlation is estimated using numpy.correlate(x,y,'full').
    Normalisation is handled by this function using the following cases:

        * 'biased': Biased estimate of the cross-correlation function
        * 'unbiased': Unbiased estimate of the cross-correlation function
        * 'coeff': Normalizes the sequence so the autocorrelations at zero
           lag is 1.0.

    :return:
        * a numpy.array containing the cross-correlation sequence (length 2*N-1)
        * lags vector

    .. note:: If x and y are not the same length, the shorter vector is
        zero-padded to the length of the longer vector.
    """
    x2 = x - np.mean(x)
    y2 = y - np.mean(y)
    N = len(x2)
    if y is None:
        y = x2
    assert len(x2) == len(y2), 'x and y must have the same length. Add zeros if needed'
    assert maxlags <= N, 'maxlags must be less than data length'

    if maxlags is None:
        maxlags = N - 1
        lags = np.arange(0, 2 * N - 1)
    else:
        assert maxlags < N
        lags = np.arange(N - maxlags - 1, N + maxlags)

    res = np.correlate(x2, y2, mode='full')

    if norm == 'biased':
        res = res[lags] / float(N)  # do not use /= !!
    elif norm == 'unbiased':
        res = res[lags] / (float(N) - abs(np.arange(-N + 1, N)))[lags]
    elif norm == 'coeff':
        Nf = float(N)
        rms = rms_flat(x2) * rms_flat(y2)
        res = res[lags] / rms / Nf
    else:
        res = res[lags]

    return res


def getTarget(data, cutCC=95, midFactor=8, start=None, end=None, mode='mean', plot=True):
    """ Returns a mean as target from data using a cutoff on the CC value.
    :param data: Images object
    :param cutCC: The percent above which to take time points
    :param midFactor: the factor around the mean to take the internal mean from
    :param start: if midFactor is 0 use start as the first index
    :param end: if midFactor is 0 use end as the last index
    :param mode: 'mean' - returns mean; 'index' returns the indices;
    :param plot: whether to plot the CC as a function of volume index
    :return: numpy array of the mean image
    """
    MeanSz = data.shape[1:]
    nRecords = data.shape[0]
    if midFactor != 0 and not (start and end):
        Mid = round(nRecords / 2)
        start = max([round(Mid - Mid / midFactor), 0])
        end = min([round(Mid + Mid / midFactor), nRecords])
    start = int(max([start, 0]))
    end = int(min([end, nRecords]))
    if midFactor == 1:
        MidData = data
    else:
        logger.info('Start' + str(start) + ' end: ' + str(end))
        if len(MeanSz) == 2:
            MidData = data[start:end, :, :]
        elif len(MeanSz) == 3:
            MidData = data[start:end, :, :, :]
    Mean = MidData.mean().toarray()

    if len(MeanSz) == 2:
        MeanVec = Mean.reshape(1, MeanSz[0] * MeanSz[1])
        CC = MidData.map(lambda vol: np.corrcoef(vol.reshape(1, MeanSz[0] * MeanSz[1]),
                                                 MeanVec)[0, 1]).toarray()
    elif len(MeanSz) == 3:
        MeanVec = Mean.reshape(1, MeanSz[0] * MeanSz[1] * MeanSz[2])
        CC = MidData.map(lambda vol: np.corrcoef(vol.reshape(1, MeanSz[0] * MeanSz[1] * MeanSz[2]),
                                                 MeanVec)[0, 1]).toarray()
    cut = np.percentile(CC, cutCC)
    Ir = (CC > cut).nonzero()[0]  # + start
    if plot:
        plt.figure()
        plt.plot(CC)
        plt.xlabel('Volume #')
        plt.ylabel('CC with mean')
        plt.title('Time points: #' + str(len(Ir)) + ', cutoff at: ' + str(cut))
    if mode == 'mean':
        if len(MeanSz) == 2:
            return MidData[Ir, :, :].mean().toarray()
        elif len(MeanSz) == 3:
            return MidData[Ir, :, :, :].mean().toarray()
    elif mode == 'index':
        return Ir
    else:
        logger.error("mode error: " + mode)


def registerByPlane(sc, data, target, upSample2=1, doShift=False, cval=0.0, oldShifts=None):
    """ for 4d array, assumes (t, x, y, z) and t could be 1.

    :param sc: SparkContext
    :param data: numpy array
    :param target: 3d target
    :param upSample2: up-sampling factor
    :param doShift: flag to preform the shift on the data
    :param cval: filling value for the shifting
    :param oldShifts: shifts from previously registration with less data
    :return: (t, 2, z) shifts
    """

    # target size
    z_num = target.shape[2]

    def volumeShiftsByPlane(target2, data2, upSample3):
        from skimage.feature import register_translation
        shifts2 = np.zeros((2, z_num))
        for i in range(z_num):
            shifts2[:, i] = register_translation(data2[:, :, i], target2[:, :, i], upSample3)[0]
        if upSample3 == 1:
            shifts2 = shifts2.astype(int)
        return shifts2

    def shiftByPlane(kv):
        from scipy.ndimage.interpolation import shift as SP_Shift
        k, v = kv
        shiftC = shiftsBC.value[int(k[0]), :, :]
        v.setflags(write=1)
        for p in range(v.shape[2]):
            v[:, :, p] = SP_Shift(v[:, :, p], shiftC[:, p], cval=cval)
        return v

    # case for one volume
    if len(data.shape) == 3 and isinstance(data, (np.ndarray, np.generic)):
        shifts = volumeShiftsByPlane(target, data, upSample2)
        if doShift:
            shiftsBC = sc.broadcast(shifts)
            return shifts, shiftByPlane((0, data))
        else:
            return shifts

    # check types and convert to thunder images
    if isinstance(data, pyspark.RDD):
        data = td.images.fromrdd(data)
    if isinstance(data, (np.ndarray, np.generic)) and sc:
        data = td.images.fromarray(data, engine=sc)
    if isinstance(data, td.base.Base):
        # check size
        if len(data.shape) == 4:
            if oldShifts is not None:
                if oldShifts.shape[0] == data.shape[0]:
                    shifts = oldShifts
                else:
                    data2 = data[oldShifts.shape[0]:, :, :, :]
                    shifts = data2.map(lambda x: volumeShiftsByPlane(x, target, upSample2), dtype=np.float32,
                                       value_shape=(2, z_num)).toarray()
                    shifts = np.concatenate((oldShifts, shifts), axis=0)
            else:
                shifts = data.map(lambda x: volumeShiftsByPlane(x, target, upSample2), dtype=np.float32,
                                  value_shape=(2, z_num)).toarray()
            if len(shifts.shape) == 2:
                shifts = shifts[:, :, np.newaxis]
            if doShift:
                shiftsBC = sc.broadcast(shifts)
                return shifts, data.map(shiftByPlane, dtype=np.int16, value_shape=data.shape[1:], with_keys=True)
            else:
                return shifts
        else:
            raise RuntimeError('shape not compatible %s' % data.shape)
    else:
        raise TypeError('type not convertible to Thunder images' + str(type(data)))


def setDefaultVal(myDict, name, value):
    """ checks if a value exists in the dictionary and returns it.
    if there isn't a value for that property sets it to value and.

    :param myDict: a dictionary we what to update
    :param name: the name of the property
    :param value: the default vale if it is not already set
    :return: the value wither default or from myDict
    """
    if name in myDict:
        return myDict[name]
    else:
        myDict[name] = value
        return value


def regroupRDDs(rdd, numGroup=10):
    """ regroup an rdd using a new key added that is 0-numGtoup-1

    :param rdd: input rdd as a (k,v) pairs
    :param numGroup: number of groups to concatenate to
    :return: a new rdd in the form of (groupNum, list of (k, v) in that group) pairs
    """
    rdd = rdd.map(lambda kv: (kv[0] % numGroup, (kv[0], kv[1])), preservesPartitioning=True)
    return rdd.groupByKey().mapValues(list)


def getTiffInfo(path):
    """ Opens the first tif file and gets the tiff header as string

    :param path: folder to tif files
    :return: volume rate in Hz
    """
    # py 2/3 comp
    first_file = glob.glob(os.path.join(path, '*.tif'))[0]
    if ScanImageTiffReader is not None and ScanImageTiffReader(first_file).metadata() != '':
        string = ScanImageTiffReader(first_file).metadata()
    else:
        tfh = tifffile.TiffFile(first_file)
        # If software key is in dict tags --> SI2016
        if 'software' in tfh.pages[0].tags:
            string = tfh.pages[0].tags['software'].value.decode('utf-8')
        else:
            string = tfh.pages[0].tags['image_description'].value.decode('utf-8')
    string = " ".join(string.split()).replace('\\', ' ')
    string = string.replace(')', '')
    string = string.replace('(', '')
    return string


def searchTiffInfo(string, search_string):
    """ Searches for data inside the tiff header string

    :param string: the tiff header string
    :param search_string: string to find in tiff header (replace space with \s)
    :return: the number after the property and the  '='
    """

    return re.search('(?<=' + search_string + '\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?', string).group(0)


def browse_data4D(data, dim1=2, dim2=3, cmap='gray', titles=None, fig_args=None, im_args=None):
    """For browsing correlation 4d matrix with shifts

    :param data: 4d numpy array
    :param dim1: the faster dimension (Z)
    :param dim2: the slower dimension (Group, Time point)
    :param cmap: the color map default to gray
    :param titles: List that will be indexed with dim2
    :param fig_args: kwargs for the plt.figure function
    :param im_args: kwargs for the plt.imshow function
    """
    s = data.shape
    otherDims = sorted(list(set(range(4)) - {dim1, dim2}))
    data = data.transpose((otherDims[0], otherDims[1], dim1, dim2))

    def view_image(index1, index2):
        if fig_args is not None:
            plt.figure(**fig_args)
        image = data[:, :, index1, index2]
        if im_args is not None:
            plt.imshow(image.T, cmap=cmap, **im_args)
        else:
            plt.imshow(image.T, cmap=cmap)
        if titles:
            plt.title(str(titles[index2]))
        plt.show()

    interact(view_image, index1=(0, s[dim1] - 1), index2=(0, s[dim2] - 1))


def browse_data4D_Dict(data, dim=2, cmap='gray', titles=None, fig_args=None, im_args=None):
    """For browsing correlation 4d matrix with shifts

    :param data: 4d numpy array
    :param dim: the faster dimension (Z)
    :param cmap: the color map default to gray
    :param titles: List that will be indexed with dim2
    :param fig_args: kwargs for the plt.figure function
    :param im_args: kwargs for the plt.imshow function
    """
    s1 = len(data)
    otherDims = sorted(list(set(range(3)) - {dim}))
    s2 = data[0].shape[dim]
    for i in range(s1):
        data[i] = data[i].transpose((otherDims[0], otherDims[1], dim))

    def view_image(index1, index2):
        if fig_args is not None:
            plt.figure(**fig_args)
        image = data[index1][:, :, index2]
        if im_args is not None:
            plt.imshow(image.T, cmap=cmap, **im_args)
        else:
            plt.imshow(image.T, cmap=cmap)
        if titles:
            plt.title(str(titles[index2]))
        plt.show()

    interact(view_image, index1=(0, s1 - 1), index2=(0, s2 - 1))


def browse_data3D(data, dim1=2, cmap='gray', titles=None, color_bar=False, fig_args=None, im_args=None):
    """For browsing correlation 3d matrix with shifts

    :param data: 3d numpy array
    :param dim1: the dimension to scroll
    :param cmap: the color map default to gray
    :param titles: List that will be indexed with dim2
    :param color_bar: flag to show a color bar
    :param fig_args: kwargs for the plt.figure function
    :param im_args: kwargs for the plt.imshow function
    """
    s = data.shape
    otherDims = sorted(list(set(range(3)) - {dim1}))
    data = data.transpose((otherDims[0], otherDims[1], dim1))

    def view_image(index1):
        if fig_args is not None:
            plt.figure(**fig_args)
        image = data[:, :, index1]
        if im_args is not None:
            plt.imshow(image.T, cmap=cmap, **im_args)
        else:
            plt.imshow(image.T, cmap=cmap)
        if color_bar:
            plt.colorbar()
        if titles:
            plt.title(str(titles[index1]))
        plt.show()

    interact(view_image, index1=(0, s[dim1] - 1))


def getBestCrop(sc, target, maxShifts, nanThreshold=200, inPaint=True, checkSize=True, origFactor=2):
    """

    :param sc: spark context
    :param target: the volume to find the best crop for
    :param maxShifts: what were the max shifts to narrow the search
    :param nanThreshold: max nan values in the cropped volume
    :param inPaint: 0=none, 1=inpaint_biharmonic, 2=nan_to_num
    :param checkSize: flag to check we don't get a crop that is smaller then current size diveded by origFactor
    :param origFactor: by how much the original volume was expanded in x and y
    :return:
    """
    sz = target.shape
    xLim = int(maxShifts[0])
    yLim = int(maxShifts[1])
    targetBC = sc.broadcast(target)
    rdd = sc.parallelize(list(product(range(yLim), range(yLim))), sc.defaultParallelism * 2)
    rdd.cache()
    rdd.count()
    borderInit = np.array([sz[0] / 4, sz[1] / 4])

    def getBestShifts(yLimits):
        bestVol = 0
        k = yLimits[0]
        m = yLimits[1]
        bestShifts = (-1, -1, k, m)
        bestNans = 0
        for i in range(xLim):
            for j in range(xLim):
                borderCurr = [borderInit[0] - i, borderInit[0] - j, borderInit[1] - k, borderInit[1] - m]
                if np.any(np.array(borderCurr) < 0):
                    break
                borderCurr = np.array(borderCurr).astype(int)
                temp = targetBC.value[borderCurr[0]:-borderCurr[1], borderCurr[2]:-borderCurr[3], :]
                curVol = temp.shape[0] * temp.shape[1]
                if curVol <= bestVol:
                    continue
                mask = np.isnan(temp)
                if np.sum(mask) > 0:
                    flag = False
                    labels = label(mask)
                    mostContigNaNs = 0
                    for w in range(1, np.max(labels) + 1):
                        temp2 = np.sum(labels == w)
                        if temp2 > nanThreshold:
                            flag = True
                            break
                        else:
                            mostContigNaNs = max(mostContigNaNs, temp2)
                    if flag:
                        break
                else:
                    mostContigNaNs = 0
                if mostContigNaNs > nanThreshold:
                    break
                if curVol > bestVol:
                    bestVol = curVol
                    bestShifts = borderCurr
                    bestNans = mostContigNaNs
        return bestShifts, bestVol, bestNans

    Flag = True
    while Flag:
        returned = rdd.map(getBestShifts).collect()
        shiftIndex = np.argmax(np.array([s[1] for s in returned]))
        borderFinal = returned[shiftIndex][0]
        if borderFinal[0] is not -1:
            Flag = False
        else:
            nanThreshold += nanThreshold * 2
            logger.info('New nanThreshold: %d' % nanThreshold)

    # make sure it is not smaller then original image
    if checkSize:
        while sz[0] - borderFinal[0] - borderFinal[1] < sz[0] / origFactor:
            borderFinal[0] -= 1
            borderFinal[1] -= 1
            logger.info(borderFinal)
        while sz[1] - borderFinal[2] - borderFinal[3] < sz[1] / origFactor:
            borderFinal[2] -= 1
            borderFinal[3] -= 1
            logger.info(borderFinal)
    borderFinal = list(map(lambda x: int(x), borderFinal))
    target = target[borderFinal[0]:-borderFinal[1], borderFinal[2]:-borderFinal[3], :]
    if inPaint == 1:
        if sz[2] > 100:
            target2 = copy.deepcopy(target).transpose(2, 0, 1)

            def inPaintPar(volume):
                mask = np.isnan(volume)
                if np.sum(mask):
                    try:
                        return inpaint_biharmonic(volume, mask)
                    except ValueError:
                        return np.nan_to_num(volume)
                else:
                    return volume

            targetImages = td.images.fromarray(target2, engine=sc).map(inPaintPar).toarray()
            target = targetImages.transpose(1, 2, 0)
        else:
            for i in range(sz[2]):
                mask = np.isnan(target[:, :, i])
                if np.sum(mask):
                    target[:, :, i] = inpaint_biharmonic(target[:, :, i], mask)
    elif inPaint == 2:
        target = np.nan_to_num(target)
    targetBC.unpersist()
    return returned, target


def add_scalebar(ax, matchx=True, matchy=True, hidex=True, hidey=True, **kwargs):
    """ Add scalebars to axes
    Adds a set of scale bars to *ax*, matching the size to the ticks of the plot
    and optionally hiding the x and y axes
    - ax : the axis to attach ticks to
    - matchx,matchy : if True, set size of scale bars to spacing between ticks
                    if False, size should be set using sizex and sizey params
    - hidex,hidey : if True, hide x-axis and y-axis of parent
    - **kwargs : additional arguments passed to AnchoredScaleBars
    Returns created scalebar object
    """

    def f(axis):
        l = axis.get_majorticklocs()
        return len(l) > 1 and (l[1] - l[0])

    if matchx:
        kwargs['sizex'] = f(ax.xaxis)
        kwargs['labelx'] = str(kwargs['sizex'])
    if matchy:
        kwargs['sizey'] = f(ax.yaxis)
        kwargs['labely'] = str(kwargs['sizey'])

    sb = AnchoredScaleBar(ax.transData, **kwargs)
    ax.add_artist(sb)

    if hidex: ax.xaxis.set_visible(False)
    if hidey: ax.yaxis.set_visible(False)

    return sb


class AnchoredScaleBar(AnchoredOffsetbox):
    def __init__(self, transform, sizex=0, sizey=0, labelx=None, labely=None, loc=4,
                 pad=0.1, borderpad=0.1, sep=2, prop=None, **kwargs):
        """
        Draw a horizontal and/or vertical  bar with the size in data coordinate
        of the give axes. A label will be drawn underneath (center-aligned).
        - transform : the coordinate frame (typically axes.transData)
        - sizex,sizey : width of x,y bar, in data units. 0 to omit
        - labelx,labely : labels for x,y bars; None to omit
        - loc : position in containing axes
        - pad, borderpad : padding, in fraction of the legend font size (or prop)
        - sep : separation between labels and bars in points.
        - **kwargs : additional arguments passed to base class constructor
        """
        from matplotlib.patches import Rectangle
        from matplotlib.offsetbox import AuxTransformBox, VPacker, HPacker, TextArea
        bars = AuxTransformBox(transform)
        if sizex:
            bars.add_artist(Rectangle((0, 0), sizex, 0, fc="none"))
        if sizey:
            bars.add_artist(Rectangle((0, 0), 0, sizey, fc="none"))

        if sizex and labelx:
            bars = VPacker(children=[bars, TextArea(labelx + 's', minimumdescent=False)],
                           align="center", pad=0, sep=sep)
        if sizey and labely:
            bars = HPacker(children=[TextArea(labely + 'DF/F'), bars],
                           align="center", pad=0, sep=sep)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=bars, prop=prop, frameon=False, **kwargs)


def getStructure(pixelSize, radius):
    """ returns a structure element to be used in filtering

    :param pixelSize: array like triplet for x, y, z size in um
    :param radius: where to reach in um
    :return: structure with 1 if dist < radius
    """
    pixelSize = np.array(pixelSize)
    w = np.round(radius / pixelSize).astype('int16')
    X, Z, Y = np.array(np.meshgrid(range(-w[0], w[0] + 1), range(-w[2], w[2] + 1), range(-w[1], w[1] + 1)))
    d = (((X * pixelSize[0]) ** 2 + (Y * pixelSize[1]) ** 2 + (Z * pixelSize[2]) ** 2) ** 0.5)
    structure = (d <= radius).astype(int)
    return structure


def showCellLocations(session):
    x = session.Sp['All']['x']
    y = session.Sp['All']['y']
    z = session.Sp['All']['z']
    for i in range(len(np.where(session.fieldMask)[0])):
        print('Cell #%d:\tx:%d, y:%d, z:%d' % (i, x[np.where(session.fieldMask)[0][i]] / 0.55,
                                               y[np.where(session.fieldMask)[0][i]] / 0.55,
                                               z[np.where(session.fieldMask)[0][i]] / 1.6))


def log_progress(sequence, every=None, size=None, name='Items'):
    """ from https://github.com/alexanderkuk/log-progress

    :param sequence: to go over and report progress
    :param every: every how many iteration to update
    :param size: size of loop if given an iterator
    :param name: to add to a label
    """
    from ipywidgets import IntProgress, HTML, VBox
    from IPython.display import display

    is_iterator = False
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = int(size / 200)  # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    display(box)

    index = 0
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{name}: {index} / ?'.format(
                        name=name,
                        index=index
                    )
                else:
                    progress.value = index
                    label.value = u'{name}: {index} / {size}'.format(
                        name=name,
                        index=index,
                        size=size
                    )
            yield record
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'success'
        progress.value = index
        label.value = "{name}: {index}".format(
            name=name,
            index=str(index or '?')
        )


def update_progress(key='1ENbQe4QtUZ6as9y8Kcr_sHl3YpwQoAs_EScIIURefSQ', gid='78168970',
                    base='/groups/svoboda/svobodalab/users/Aaron',
                    json_keyfile='/groups/svoboda/svobodalab/users/moharb/python-aed5a054cd02.json'):
    """ see https://github.com/burnash/gspread

    :param key: sheet ket
    :param gid: sheet gid
    :param base: base folder for raw data
    :param json_keyfile: for write access
    :return:
    """
    # google sheets to pandas dataframe
    response = requests.get('https://docs.google.com/spreadsheet/ccc?key=' + key + '&output=csv&gid=' + gid)
    assert response.status_code == 200, 'Wrong status code'
    f = io.StringIO(response.content.decode('utf-8'))
    sessionsDF = pd.read_csv(f)

    # clean up dataframe
    wr_not_nan = pd.notnull(sessionsDF.WR)
    m_roi = sessionsDF.Type == 'MROI'
    sessions_clean = sessionsDF[wr_not_nan & m_roi]
    # del sessions_clean['Animal']
    # del sessions_clean['Type']
    # del sessions_clean['StartWeight']
    # del sessions_clean['EndWeight']
    # del sessions_clean['Stack']
    # del sessions_clean['Pitch']
    # del sessions_clean['Roll']

    # OAth to write back
    scope = ['https://spreadsheets.google.com/feeds']
    credentials = ServiceAccountCredentials.from_json_keyfile_name(json_keyfile, scope)
    gc = gspread.authorize(credentials)
    sh = gc.open_by_key(key)
    worksheet = sh.worksheet("Imaging")
    indexes = sessions_clean.index.values
    for index in log_progress(indexes):
        f = sessions_clean.loc[index]
        date = pd.to_datetime(f['StartTime']).strftime('%y%m%d')
        run = 'Run' + str(int(f['Run']))
        animal = f['WR']
        base_directory = os.path.join(base, '%s_%s' % (date, animal))
        fov = 'FOV' + str(int(f['FOV']))
        full_directory = os.path.join(base_directory, run)
        database_directory = os.path.join(base, 'Database', animal, fov, date + run)
        row = index + 2

        # Select a range
        cell_list = worksheet.range('S%d:Z%d' % (row, row))
        values = list()
        # raw data in dm11
        values.append(os.path.isdir(full_directory))
        # Sp.mat copied over
        values.append(os.path.isfile(os.path.join(full_directory, 'Sp.mat')))
        # done with registration
        values.append(os.path.isfile(os.path.join(base_directory, run + 'step5.p')))
        # copied over to Database
        values.append(os.path.isfile(os.path.join(database_directory, 'Sp.mat')) and os.path.isfile(
            os.path.join(database_directory, 'expended_new.tif')))
        # prepared to trace
        values.append(os.path.isfile(os.path.join(database_directory, 'Sp_2.mat')) and os.path.isfile(
            os.path.join(database_directory, 'prepareMasksAuto.mat')) and os.path.isfile(
            os.path.join(database_directory, 'prepareMasksAutoSWC.mat')) and os.path.isfile(
            os.path.join(database_directory, 'Session.tif')) and os.path.isfile(
            os.path.join(database_directory, 'top_vol_up.mat')))
        # traced dendrite
        values.append(os.path.isfile(os.path.join(database_directory, 'dendrite.swc')))
        # traced spines
        values.append(os.path.isfile(os.path.join(base_directory, run + 'mask.mat')))
        # made time2 new
        values.append(os.path.isfile(os.path.join(base_directory, run + 'time2_new.p')))
        # update values
        for cell, value in zip(cell_list, values):
            cell.value = value
        # Update in batch
        worksheet.update_cells(cell_list)


def shifts_projection(sc, clean):
    """ Register a volume in 3d by projecting and using CC on the projections to get shifts
    Averages the x (y, z) shifts for the xy and xz projection results

    :param sc: Spark context
    :param clean: 3d input data
    :return:  registratred images object
    """
    def shifts_projected(clean, axis):
        projected = clean.map(lambda x: x.mean(axis=axis)[:, :, np.newaxis])
        target = getTarget(projected, 30, 1)
        shifts = registerByPlane(sc, projected, target[:, :, np.newaxis], 10, False)
        return shifts[:, :, 0]

    # shifts_xy = shifts_projected(clean, 2)
    shifts_xz = shifts_projected(clean, 1)
    shifts_yz = shifts_projected(clean, 0)

    # x_shifts = np.mean(np.stack((shifts_xz[:, 0], shifts_xy[:, 0])), axis=0)
    z_shifts = np.mean(np.stack((shifts_xz[:, 1], shifts_yz[:, 1])), axis=0)
    # y_shifts = np.mean(np.stack((shifts_yz[:, 0], shifts_xy[:, 1])), axis=0)
    plt.figure()
    plt.plot(shifts_xz[:, 1])
    plt.plot(shifts_yz[:, 1])
    plt.plot(z_shifts)
    plt.title('Z')
    # plt.figure()
    # plt.plot(shifts_xz[:, 0])
    # plt.plot(shifts_xy[:, 0])
    # plt.plot(x_shifts)
    # plt.title('X')
    # plt.figure()
    # plt.plot(shifts_yz[:, 0])
    # plt.plot(shifts_xy[:, 1])
    # plt.plot(y_shifts)
    # plt.title('Y')
    # shifts_all = np.stack((x_shifts, y_shifts, z_shifts))

    def initReg(kv):
        from scipy.ndimage.interpolation import shift
        index, volume = kv
        current_shift = (0, 0, -1 * z_shifts[int(index[0])])
        shifted = shift(volume, current_shift)
        return shifted.astype(np.int16)

    reg = clean.map(initReg, with_keys=True, value_shape=clean.shape[1:], dtype=np.int16)
    reg.cache()
    reg.count()
    return reg


def count_timepoints(sc, session, files):
    """ counts timepoints in each tif file, tries to use ScanImage reader for speed, falls back to tifffile
        assumes the first dimension is number of time points in each file
    :param sc: Spark Context
    :param session: Session object
    :param files: list of files
    :return:
    """
    tuples = zip(range(len(files)), files)
    files_sc = sc.parallelize(tuples)

    def count_planes(kv):
        index, path2 = kv
        try:
            from ScanImageTiffReader import ScanImageTiffReader
            img = ScanImageTiffReader(path2).data()
        except Exception:
            import tifffile
            img = tifffile.imread(path2)
        return img.shape[0]

    data2 = files_sc.map(count_planes).collect()
    frame_numbers = np.array(data2)
    vol_numbers = frame_numbers / len(session.fieldMask)
    return vol_numbers.astype(int)


def separate_sessions(sc, session, reg_data, run_names=('Run1', 'Run2'), run_number=(1, 2)):
    """ Two runs that were registered and time course (ed) together could be separated

    :param sc: spark context
    :param session: the session to split
    :param reg_data: registered images object to split
    :param run_names: the names of the runs as a list example: ['Run1','Run2']
    :param run_number: a list of the numbers (1, 2)
    :return:
    """
    from prep.VisSession import VisSession
    from prep.IO import saveRegData
    # separate
    start = 0
    for run, name in zip(run_number, run_names):
        files = glob.glob(os.path.join(session.path, '') + name + '*.tif')
        timepoints = count_timepoints(sc, session, files)
        stop = start + sum(timepoints)
        logger.info('Run %s, Number %d, start: %d, stop %d' % (name, run, start, stop))
        session_temp = VisSession(animalID=session.animalID, date=session.date, run='Run' + str(run))
        session_temp.initBase(sc)
        session_temp.start = 0
        session_temp.stop = len(files)
        reg_temp = reg_data[start:stop, :, :, :]

        saveRegData(session_temp, reg_temp, overwrite=False)
        reg_dict_temp = copy.deepcopy(session.regDict)
        for key, value in iteritems(session.regDict):
            if isinstance(value, np.ndarray):
                shape = value.shape
                loc = np.where(np.array(shape) == reg_data.shape[0])[0]
                if len(loc) > 0:
                    print('%s: %s - %s' % (key, shape, loc[0]))
                    if loc[0] == 0:
                        reg_dict_temp[key] = value[start:stop, ...]
        session_temp.regDict = reg_dict_temp
        time_dict_temp = copy.deepcopy(session.timeDict)
        for key, value in iteritems(session.timeDict):
            if isinstance(value, np.ndarray):
                shape = value.shape
                loc = np.where(np.array(shape) == reg_data.shape[0])[0]
                if len(loc) > 0:
                    print('%s: %s - %s' % (key, shape, loc[0]))
                    if loc[0] == 0:
                        time_dict_temp[key] = value[start:stop, ...]
                    if loc[0] == 1:
                        time_dict_temp[key] = value[:, start:stop, ...]
        session_temp.timeDict = time_dict_temp
        session_temp.save('time2_new')
        start = stop


def rotate(offset, points, angle):
    theta = np.radians(angle)
    c, s = np.cos(theta), np.sin(theta)
    r = np.matrix([[c, -s], [s, c]])
    return np.array(np.dot(points - offset, r) + offset)


def downsample(array, factor=2):
    """ from https://stackoverflow.com/questions/30379311/fast-way-to-take-average-of-every-n-rows-in-a-npy-array

    :param array: numpy array
    :param factor: factor to donw sample by averaging
    :return:
    """
    cum = np.nancumsum(array, 0)
    result = cum[factor - 1::factor] / float(factor)
    result[1:] = result[1:] - result[:-1]

    remainder = array.shape[0] % factor
    if remainder != 0:
        if remainder < array.shape[0]:
            lastAvg = (cum[-1]-cum[-1-remainder])/float(remainder)
        else:
            lastAvg = cum[-1]/float(remainder)
        result = np.vstack([result, lastAvg])
    return result


def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))


def length(v):
    return math.sqrt(dotproduct(v, v))


def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
