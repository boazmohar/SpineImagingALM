# Encoding: utf-8
""" Steps module for processing Thunder Images object with traceability and repeatability
 steps follow Step Class interface conventions
"""
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import fourier_shift, zoom
from scipy.ndimage.morphology import white_tophat
from scipy.stats.mstats import zscore
from skimage.feature import register_translation

from prep.Utils import setDefaultVal


def filterDataStep(**kwargs):
    """ filter data time point

    :param kwargs: either provide a start and stop or an index of wanted time points
    :return: cropped data
    """
    data = kwargs['data']
    start = setDefaultVal(kwargs, 'start', None)
    stop = setDefaultVal(kwargs, 'stop', None)
    index = setDefaultVal(kwargs, 'index', None)
    if index:
        return data[index, :, :, :]
    elif 0 <= start < stop:
        return data[slice(start, stop), :, :, :]
    else:
        print('error in input to filterDataStep')


def cleanDataStep(**kwargs):
    """    Returns an Images object after zeroing all pixel values lower then cutoff
    data: Images object
    cutoff: Pixel value underneath all values will be zero
    :return: a new Image object
    """
    cutoff = kwargs['cutoff']
    data = kwargs['data']

    def remNoise(x):
        x[x < cutoff] = 0
        return x

    return data.map(remNoise)


def zoomDataStep(**kwargs):
    """    Returns an Images object after zeroing all pixel values lower then cutoff
    data: Images object
    cutoff: Pixel value underneath all values will be zero
    :return: a new Image object
    """
    zoom_val = kwargs['zoom']
    data = kwargs['data']

    return data.map(lambda x: zoom(x, zoom_val))


def cropDataStep(**kwargs):
    """    Returns an Images object after cropping in X, Y and Z
    data: Images object
    param: xStart,xStop: indexes to start and end pixel cropping
    param: yStart,yStop indexes to start and end line cropping
    both default to: start --> 0, end --> size X/Y
    param: planes: index of planes in Z to take

    :return: a new Image object
    """
    data = kwargs['data']
    sz = data.shape

    tIndex = setDefaultVal(kwargs, 'tIndex', None)
    xStart = setDefaultVal(kwargs, 'xStart', 0)
    yStart = setDefaultVal(kwargs, 'yStart', 0)
    xStop = setDefaultVal(kwargs, 'xStop', sz[2])
    yStop = setDefaultVal(kwargs, 'yStop', sz[1])
    planes = setDefaultVal(kwargs, 'planes', None)
    if tIndex is not None:
        data = data[np.where(tIndex)[0], :, :, :]
    if planes is not None:
        return data[:, slice(yStart, yStop), slice(xStart, xStop), np.where(planes)[0]]
    else:
        return data[:, slice(yStart, yStop), slice(xStart, xStop), :]


def cropDictStep(**kwargs):
    """

    Parameters
    ----------
    kwargs

    Returns
    -------

    """
    data = kwargs['data']
    sz = data.shape
    cropDict = kwargs['cropDict']
    tIndex = setDefaultVal(cropDict, 'tIndex', None)
    xStart = setDefaultVal(cropDict, 'xStart', 0)
    yStart = setDefaultVal(cropDict, 'yStart', 0)
    xStop = setDefaultVal(cropDict, 'xStop', sz[2])
    yStop = setDefaultVal(cropDict, 'yStop', sz[1])
    zIndex = setDefaultVal(cropDict, 'zIndex', None)
    if tIndex is not None:
        data = data[tIndex, :, :, :]
    if zIndex is not None:
        data = data[:, :, :, np.where(zIndex)[0]]
    return data[:, slice(yStart, yStop), slice(xStart, xStop), :]


def scanPhaseStep(**kwargs):
    """corrects for bi-directional scan phase errors by computing the median of median shifts
    across the data set and applying it to all the data

    Parameters
    ----------
    data: Images RDD with bidirectional scan phase problem
    subPixel: (default: 50) up-sampled matrix-multiplication DFT
                     to achieve arbitrary subpixel precision will 20 will
                     have 1/20 pixel precision
    timepoints: (default: 5000) how many time points to use

    Returns
    -------
    corrected Images RDD
    """

    data = kwargs['data']
    subPixel = setDefaultVal(kwargs, 'subPixel', 50)
    threshold = setDefaultVal(kwargs, 'threshold', np.nan)
    timepoints = setDefaultVal(kwargs, 'timepoints', 5000)
    fixPerTimepoint = setDefaultVal(kwargs, 'fixPerTimepoint', False)
    print(threshold)
    def medianShift(image):
        newShift = np.zeros(Sz[3])
        for i in range(Sz[3]):
            evenLines = image[range(0, Sz[1], 2), :, i]
            oddLines = image[range(1, Sz[1], 2), :, i]
            minLines = min([evenLines.shape[0], oddLines.shape[0]])
            cropEven = evenLines[:minLines, :]
            cropOdd = oddLines[:minLines, :]
            shift = register_translation(cropEven, cropOdd, subPixel)
            newShift[i] = shift[0][1]
        if not np.isnan(threshold):
            Means = image.mean(axis=(0, 1))
            return np.median(newShift[Means > threshold])
        else:
            return np.median(newShift)

    def fixImage(image):
        for i in range(Sz[3]):
            oddLines = image[range(1, Sz[1], 2), :, i]
            Shift = (0, medianShift)
            shiftedOdd = fourier_shift(np.fft.fftn(oddLines), Shift)
            shiftedOdd = np.real(np.fft.ifftn(shiftedOdd))
            image[range(1, Sz[1], 2), :, i] = shiftedOdd
        return image

    def fixPerTimepointFunc(kv):
        k, image = kv
        shift = (0, shifts[k])
        for i in range(Sz[3]):
            oddLines = image[range(1, Sz[1], 2), :, i]
            shiftedOdd = fourier_shift(np.fft.fftn(oddLines), shift)
            shiftedOdd = np.real(np.fft.ifftn(shiftedOdd))
            image[range(1, Sz[1], 2), :, i] = shiftedOdd
        return image

    Sz = data.shape
    if Sz[0] > timepoints and not fixPerTimepoint:
        rand = np.sort(np.random.permutation(range(Sz[0]))[:timepoints])
        data2 = data[rand, :, :, :]
    else:
        data2 = data
    shifts = data2.map(lambda x: medianShift(x)).toarray()
    medianShift = np.median(shifts)
    plt.hist(shifts, np.max([np.ceil(Sz[0] / 1000.0), 10]).astype(int))
    plt.title('shifts median: %0.2f px' % medianShift)
    if fixPerTimepoint:
        return data.map(lambda x: fixPerTimepointFunc(x), with_keys=True)
    else:
        return data.map(lambda x: fixImage(x))


def zScoreStep(**kwargs):
    """

    Parameters
    ----------
    kwargs: axis: on which axis to apply the z-score, default to None

    Returns
    -------
    data after z-score
    """
    data = kwargs['data']
    axis = setDefaultVal(kwargs, 'axis', None)
    return data.map(lambda x: zscore(x, axis=axis))


def topHatStep(**kwargs):
    """

    Parameters
    ----------
    kwargs: structureSize: size in pixels of the structure to preform the top-hat filter

    Returns
    -------
    data after top hat filtering of each plane in the volume
    """

    data = kwargs['data']
    structureSize = setDefaultVal(kwargs, 'structureSize', 3)

    def topHatVolume(x):
        for i in range(x.shape[2]):
            x[:, :, i] = white_tophat(x[:, :, i], structureSize)
        return x
    return data.map(topHatVolume)
