# Encoding: utf-8
""" Embedding module to calculate embedding transform, get example data and browse it

"""
from __future__ import print_function

import copy
import logging
import sys
import time
from builtins import map

import matplotlib.pyplot as plt
import numpy as np
import pyspark
import scipy.signal as signal
import thunder as td

# Setup logging
logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def initEmbedDict(session, downSample=1, zStep=1.0, ZExpandPerUm=0.0, useRealZ=True, ZLPfilter=4, zFWHM=2.5):
    """ initialize b dictionary for shift and corr registration by volume

    :param session: current session object
    :param downSample: down sample factor for dense stack
    :param zStep: step for Z in um
    :param ZExpandPerUm: expansion factor to correct for remote focusing (old value was 0.0004375 04272016)
    :param useRealZ: if True will correct each line for error recorded from Z piezo/voicecoil
    :param ZLPfilter: number of planes to use Z filtering
    :param zFWHM: in um
    :return: registration dictionary
    """
    # setting properties from input vars
    embedDict = dict()
    otherVars = locals()
    del otherVars['session']
    embedDict.update(otherVars)

    # transferring needed data from session
    embedDict['xSizeOrig'] = session.xSizeOrig
    embedDict['ySizeOrig'] = session.ySizeOrig
    embedDict['xSize'] = session.xSize / downSample
    embedDict['ySize'] = session.ySize / downSample
    embedDict['UM_1X'] = session.UM_1X
    embedDict['path'] = session.path
    embedDict['ZActual'] = session.ZActual
    embedDict['xScale'] = session.xScale
    embedDict['yScale'] = session.yScale
    embedDict['xScaleAnat'] = session.xScaleAnat
    embedDict['yScaleAnat'] = session.yScaleAnat
    embedDict['x'] = session.x
    embedDict['y'] = session.y
    embedDict['z'] = session.z
    embedDict['TrueFieldsCenterSamp'] = session.TrueFieldsCenterSamp

    return embedDict


def prepareCoordinates(embedDict, subMin=True):
    """ converts the x,y,z coordinates of the volume to brain space
    if embedDict['useRealZ'] considers the error from measuring the Z feedback
    :param embedDict: the embedding dictionary
    :param subMin: subtract the min value from each dimension.
    :return: brainFields - 3 by #pixels array with x,y,z discrete location in brain space

    """
    x = embedDict['x']
    y = embedDict['y']
    z = embedDict['z']
    zStep = embedDict['zStep']
    xSizeOrig = embedDict['xSizeOrig']
    ySizeOrig = embedDict['ySizeOrig']
    xScaleAnat = embedDict['xScaleAnat']
    yScaleAnat = embedDict['yScaleAnat']
    ImagingPixels = embedDict['xSize']
    ImagingLines = embedDict['ySize']
    ZActual = embedDict['ZActual']
    TrueFieldsCenterSamp = embedDict['TrueFieldsCenterSamp']
    UM_1X = embedDict['UM_1X']
    ZExpandPerUm = embedDict['ZExpandPerUm']
    nFields = len(x)
    xSize = float(xSizeOrig) * float(embedDict['xScale'])
    ySize = float(ySizeOrig) * float(embedDict['yScale'])

    xStep = float(xSize) / float(ImagingPixels)
    yStep = float(ySize) / float(ImagingLines)

    centIdx = int(round(ImagingLines / float(2)))
    ZIdx = np.asarray(range(0, int(ImagingLines)), dtype='int64') - centIdx
    # convert to line-by-line positions
    zLines = np.zeros((int(ImagingLines), nFields))
    xLines = np.zeros((int(ImagingLines), nFields))
    yLines = np.zeros((int(ImagingLines), nFields))

    N = 2  # Filter order
    Wn = float(1) / (ImagingLines / 2)  # Cutoff frequency at half field
    B, A = signal.butter(N, Wn, output='ba')

    if embedDict['useRealZ']:
        for i in range(0, nFields):
            # smooth
            Zfilter = copy.deepcopy(ZActual[(TrueFieldsCenterSamp[i] + ZIdx)])
            Zfilter = signal.filtfilt(B, A, Zfilter)
            zLines[:, i] = Zfilter
            zScaleFactor = (float(1) + (zLines[:, i] * ZExpandPerUm))
            zScaleOffset = ((zScaleFactor * UM_1X) - UM_1X) / 2
            xLines[:, i] = x[i] * xScaleAnat * zScaleFactor - (xSize / float(2)) * zScaleFactor + zScaleOffset
            yLines[:, i] = y[i] * yScaleAnat * zScaleFactor - (ySizeOrig / float(2)) * zScaleFactor + np.asarray(
                range(0, int(ImagingLines)), dtype='float32') * yStep + zScaleOffset
    else:
        z = np.int16(np.round(z / zStep))
        for i in range(0, nFields):
            zLines[:, i] = z[i]
            zScaleFactor = (float(1) + (zLines[:, i] * ZExpandPerUm))
            zScaleOffset = ((zScaleFactor * UM_1X) - UM_1X) / 2
            xLines[:, i] = x[i] * xScaleAnat * zScaleFactor - (xSize / float(2)) * zScaleFactor + zScaleOffset
            yLines[:, i] = y[i] * yScaleAnat * zScaleFactor - (ySizeOrig / float(2)) * zScaleFactor + np.asarray(
                range(0, int(ImagingLines)), dtype='float32') * yStep + zScaleOffset

    if subMin:
        zLines -= np.min(zLines)
        xLines -= np.min(xLines)
        yLines -= np.min(yLines)

    zLines = np.round(zLines / zStep).astype('int64')
    xLines = np.round(xLines / xStep).astype('int64')
    yLines = np.round(yLines / yStep).astype('int64')

    ImagingPixels = int(ImagingPixels)
    brainFields = np.zeros((int(ImagingLines), int(ImagingPixels), 3, nFields))
    xRange = np.asarray(range(0, int(ImagingPixels)), dtype='float32')
    for i in range(0, nFields):
        for j in range(0, int(ImagingLines)):
            brainFields[j, :, 0, i] = xLines[j, i] + np.round(xRange).astype('int64')
            brainFields[j, :, 1, i] = yLines[j, i] * np.ones((1, ImagingPixels), dtype='int64')
            brainFields[j, :, 2, i] = zLines[j, i] * np.ones((1, ImagingPixels), dtype='int64')

    embedDict['brainFields'] = brainFields


def prepareCoordinatesZExp(embedDict, grpZPos, flip=False, subMin=True):
    """ same as prepareCoordinates but prepares b space in which we can embed groups that are distributed in Z
    according to grpZPos
    :param embedDict: the embedding dictionary
    :param grpZPos: positions of the groups in z (um)
    :param flip: flip the order of the groups (it is arbitrary)
    :param subMin: subtract the min value from each dimension.
    :return: brainFieldsZExp - 3 by #pixels array with x,y,z discrete location in brain space
    """

    x = embedDict['x']
    y = embedDict['y']
    z = embedDict['z']
    zStep = embedDict['zStep']
    xSizeOrig = embedDict['xSizeOrig']
    ySizeOrig = embedDict['ySizeOrig']
    xScaleAnat = embedDict['xScaleAnat']
    yScaleAnat = embedDict['yScaleAnat']
    ImagingPixels = embedDict['xSize']
    ImagingLines = embedDict['ySize']
    ZActual = embedDict['ZActual']
    ZLPfilter = embedDict['ZLPfilter']
    TrueFieldsCenterSamp = embedDict['TrueFieldsCenterSamp']
    UM_1X = embedDict['UM_1X']
    ZExpandPerUm = embedDict['ZExpandPerUm']
    nFields = len(x)
    nZShifts = grpZPos.shape[0]
    embedDict['nZShifts'] = nZShifts
    xSize = float(xSizeOrig) * float(embedDict['xScale'])
    ySize = float(ySizeOrig) * float(embedDict['yScale'])

    xStep = float(xSize) / float(ImagingPixels)
    yStep = float(ySize) / float(ImagingLines)

    centIdx = int(round(ImagingLines / float(2)))
    ZIdx = np.asarray(range(int(ImagingLines)), dtype='int64') - centIdx
    # convert to line-by-line positions
    zLines = np.zeros((int(ImagingLines), nFields * nZShifts))
    xLines = np.zeros((int(ImagingLines), nFields * nZShifts))
    yLines = np.zeros((int(ImagingLines), nFields * nZShifts))

    # Design the Buterworth filter
    N = 2  # Filter order
    Wn = float(1) / ZLPfilter  # Cutoff frequency
    B, A = signal.butter(N, Wn, output='ba')
    ZActual = signal.filtfilt(B, A, ZActual)

    if embedDict['useRealZ']:
        for i in range(0, nFields):
            for j in range(0, nZShifts):
                newFieldID = i * nZShifts + j
                zLines[:, newFieldID] = ZActual[(TrueFieldsCenterSamp[i] + ZIdx)]
                zScaleFactor = (float(1) + (zLines[:, newFieldID] * ZExpandPerUm))
                zScaleOffset = ((zScaleFactor * UM_1X) - UM_1X) / 2
                xLines[:, newFieldID] = x[i] * xScaleAnat * zScaleFactor - \
                                        (xSize / float(2)) * zScaleFactor + zScaleOffset
                yLines[:, newFieldID] = np.median(y[i] * yScaleAnat * zScaleFactor -
                                                  (ySizeOrig / float(2)) * zScaleFactor + zScaleOffset) + \
                                        np.asarray(range(int(ImagingLines)), dtype='float32') * yStep
                if flip:
                    zLines[:, newFieldID] = zLines[:, newFieldID] + grpZPos[j]
                else:
                    zLines[:, newFieldID] = zLines[:, newFieldID] - grpZPos[j]
    else:
        z = np.int16(np.round(z / zStep))
        for i in range(0, nFields):
            for j in range(0, nZShifts):
                newFieldID = i * nZShifts + j
                zLines[:, newFieldID] = z[i]
                zScaleFactor = (float(1) + (zLines[:, newFieldID] * ZExpandPerUm))
                zScaleOffset = ((zScaleFactor * UM_1X) - UM_1X) / 2
                xLines[:, newFieldID] = x[i] * xScaleAnat * zScaleFactor - (xSize / float(2)) \
                                                                           * zScaleFactor + zScaleOffset
                yLines[:, newFieldID] = np.median(y[i] * yScaleAnat * zScaleFactor - (ySizeOrig / float(2))
                                                  * zScaleFactor + zScaleOffset) + \
                                        np.asarray(range(int(ImagingLines)), dtype='float32') * yStep
                if flip:
                    zLines[:, newFieldID] = zLines[:, newFieldID] + grpZPos[j]
                else:
                    zLines[:, newFieldID] = zLines[:, newFieldID] - grpZPos[j]

    if subMin:
        zLines -= np.min(zLines)
        xLines -= np.min(xLines)
        yLines -= np.min(yLines)

    zLines = np.round(zLines / zStep).astype('int64')
    xLines = np.round(xLines / xStep).astype('int64')
    yLines = np.round(yLines / yStep).astype('int64')

    nFields = zLines.shape[1]
    brainFields = np.zeros((int(ImagingLines), int(ImagingPixels), 3, nFields))
    xRange = np.asarray(range(int(ImagingPixels)), dtype='float32')
    for i in range(nFields):
        for j in range(int(ImagingLines)):
            brainFields[j, :, 0, i] = xLines[j, i] + np.round(xRange).astype('int64')
            brainFields[j, :, 1, i] = yLines[j, i] * np.ones((1, int(ImagingPixels)), dtype='int64')
            brainFields[j, :, 2, i] = zLines[j, i] * np.ones((1, int(ImagingPixels)), dtype='int64')

    embedDict['brainFieldsZExp'] = brainFields


def getFieldTFormPar(sc, embedDict, zExpand=False, xDim=None, yDim=None, zIndex=None):
    """ Calculates the transformation from discrete brain space to expanded space

    :param sc: SparkContext
    :param embedDict: the embedding dictionary
    :param zExpand: using the Z expended transformation or not
    :param xDim: (start, stop) tuple
    :param yDim: (start, stop) tuple
    :param zIndex: list of indexes to exclude
    :return: Transformation dictionary with the number of pixels and coorArray
    """
    Start = time.time()
    zStep = embedDict['zStep']
    zFWHM = embedDict['zFWHM']
    # get cropping indices
    if zExpand:
        brainFields = embedDict['brainFieldsZExp']
    else:
        brainFields = embedDict['brainFields']
    dims = brainFields.shape
    # x dim to index
    if xDim:
        xIndex = np.zeros(dims[1], dtype=bool)
        xIndex[range(xDim[0], xDim[1])] = True
    else:
        xIndex = np.ones(dims[1], dtype=bool)
    # y dim to index
    if yDim:
        yIndex = np.zeros(dims[0], dtype=bool)
        yIndex[range(yDim[0], yDim[1])] = True
    else:
        yIndex = np.ones(dims[0], dtype=bool)
    # z index
    if zIndex is not None:
        # excluding indexes that are in zIndex
        if zExpand:
            raise NotImplementedError('index exclusion in expended not implemented')
            # zIndexSize = zIndex.shape[0]
            # nZShifts = embedDict['nZShifts']
            # for z in range(1, nZShifts):
            #     zIndex = np.concatenate((zIndex, zIndex[:zIndexSize]+nZShifts*z))
        # both for expended and not
        temp = np.ones(dims[3], dtype=bool)
        temp[zIndex] = False
        zIndex = temp
    else:
        # Using all fields
        zIndex = np.ones(dims[3], dtype=bool)
    # get cropped linearFields
    brainFields = brainFields[yIndex, :, :, :]
    brainFields = brainFields[:, xIndex, :, :]
    brainFields = brainFields[:, :, :, zIndex]
    linearFields = brainFields.transpose([0, 1, 3, 2])
    si = linearFields.shape
    linearFields = linearFields.reshape((si[0] * si[1] * si[2], si[3]), order="F").astype(np.int32)
    count = si[0] * si[1] * si[2]

    # save the results in b dict
    TFormDict = dict()
    TFormDict['coorArray'] = linearFields
    TFormDict['brainFields'] = brainFields
    TFormDict['pixelsPerVolume'] = count
    TFormDict['xIndex'] = xIndex
    TFormDict['yIndex'] = yIndex
    TFormDict['zIndex'] = zIndex

    # prepare for calculation
    embedSize = np.ceil(np.amax(linearFields, 0)) + 1
    a = linearFields[:, 0:2]
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, ia, ic = np.unique(b, return_index=True, return_inverse=True)
    allUniqueXY = a[ia]
    nUniqueXY = allUniqueXY.shape[0]
    linearFieldsBC = sc.broadcast(linearFields)

    def calcTform(key):
        """  parallel by z field transformation calculation
        :param key: current z position
        :return: Transformation for that field
        """
        from scipy.stats import norm
        key = int(key[0])
        fieldsTformInt = [None] * linearFieldsBC.value.shape[0]
        embedSize2 = embedSize.astype('int')
        for j in range(nUniqueXY):
            sameXYIdx = np.nonzero(ic == j)[0]
            zDist = np.absolute((key - 0.5) - linearFieldsBC.value[sameXYIdx, 2])
            validXYIdx = sameXYIdx[((zDist * zStep) < zFWHM).nonzero()[0]]
            validXYIdx = validXYIdx.astype('int64')
            zDist = np.absolute((key - 0.5) - linearFieldsBC.value[validXYIdx, 2])
            zDistNorm = norm.pdf((zDist * zStep) / (zFWHM / 2.355 / 2)) / norm.pdf(0)
            fracContribute = zDistNorm / zDistNorm.shape[0]
            for k in range(len(validXYIdx)):
                if not np.any(fieldsTformInt[validXYIdx[k]]):
                    fieldsTformInt[validXYIdx[k]] = np.asarray([np.ravel_multi_index(
                        [allUniqueXY[j][0].astype('int64'), allUniqueXY[j][1].astype('int64'), key], embedSize2,
                        order='F'), fracContribute[k]]).reshape(-1, 1)
                else:
                    fieldsTformInt[validXYIdx[k]] = np.concatenate(
                        (fieldsTformInt[validXYIdx[k]], np.asarray([np.ravel_multi_index(
                            [allUniqueXY[j][0].astype('int64'), allUniqueXY[j][1].astype('int64'), key],
                            embedSize2,
                            order='F'), fracContribute[k]]).reshape(-1, 1)), axis=1)
        return fieldsTformInt

    idxList = zip(range(int(embedSize[2])))
    allFieldsTform = sc.parallelize(idxList).flatMap(lambda x: calcTform(x)).collect()
    allFieldsTform2 = np.asarray(allFieldsTform)
    allFieldsTform2 = allFieldsTform2.reshape(-1, int(embedSize[2]), order='F')
    fieldsTform = [None] * allFieldsTform2.shape[0]
    for i in range(allFieldsTform2.shape[0]):
        a = allFieldsTform2[i, :]
        a = a[np.array(list(map(lambda b: type(b) == np.ndarray, a)))]
        if a.shape[0] > 0:
            fieldsTform[i] = np.concatenate(a, axis=1)
    TFormDict['embedSize'] = embedSize.astype(int)
    TFormDict['fieldsTform'] = fieldsTform
    End = time.time()
    logger.info('fieldsTform:: time(s): ' + str(End - Start))
    return TFormDict


def getExampleVol(sc, data, TFormDict, start=None, end=None, project=False, verbose=False, max_project=False):
    """ Embeds data

    :param sc: SparkContext
    :param data: Images object, equivalent rdd, or 4D numpy array (1st dimension is number of volumes)
    :param TFormDict: dictionary to take the transformation from
    :param start: Index to start from (if None will take 0)
    :param end: Index to end (if None will take nrecords)
    :param project: if True will return b Z projected embedding
    :param verbose: if True will print
    :return: 3D/4D embedded numpy array
    """
    tStart = time.time()
    if start and end and start == end:
        logger.error('start and end should be different')
        return
    if isinstance(data, pyspark.RDD):
        data = td.images.fromrdd(data)
    elif isinstance(data, (np.ndarray, np.generic)):
        # if only one time point add b dimension so Thunder would see it as 1 time point
        if len(data.shape) == 3:
            data = np.expand_dims(data, 0)
        data = td.images.fromarray(data, engine=sc)
    elif isinstance(data, td.base.Base):
        pass
    else:
        logger.error('wrong data type', type(data))
        return

    nrecords = data.shape[0]
    if not start:
        start = 0
    if not end:
        end = nrecords + 1

    fieldsTformBC = sc.broadcast(TFormDict['fieldsTform'])
    embedSize = TFormDict['embedSize'].astype(int)
    def embedVol(fields):
        """

        :param fields:
        :return:
        """

        fieldsTformInt = fieldsTformBC.value
        fields = fields.flatten(order='F')
        embeddedVol = np.empty(np.prod(embedSize))
        embeddedVol.fill(np.NAN)
        for i in range(len(fieldsTformInt)):
            Tform = fieldsTformInt[i].transpose()
            for j in range(Tform.shape[0]):
                embeddedVol[int(Tform[j, 0])] = np.nansum([embeddedVol[int(Tform[j, 0])], fields[i] * Tform[j, 1]])
        embeddedVol = embeddedVol.reshape(embedSize, order='F')
        if project:
            if max_project:
                embeddedVol = np.nanmax(embeddedVol, axis=2)
            else:
                embeddedVol = np.nanmean(embeddedVol, axis=2)
        return embeddedVol.astype(dtype='float32')

    SampleData = data[start:end]
    SampleData.cache()
    embeddedData = SampleData.map(embedVol).toarray()
    SampleData.uncache()

    tEnd = time.time()
    if verbose:
        logger.info('getExampleVol:: time(s): ' + str(tEnd - tStart))
    return embeddedData


def getExtendedVol(sc, groupZ, embedDict, TForm=None):
    """ embedding the Z groups

    :param sc: SparkContext
    :param embedDict: embedDict
    :param groupZ: the Z aggregated groups
    :param TForm: the expended transformation dict in None will calculate
    :return:
    """

    if TForm is None:
        TForm = getFieldTFormPar(sc, embedDict, zExpand=True)
    groupZ = groupZ.transpose([1, 2, 3, 0])
    si = groupZ.shape
    groupZ = groupZ.reshape(1, si[0], si[1], si[2] * si[3])
    return getExampleVol(sc, groupZ, TForm)


def getTFrom(sc, session, return_t=False):
    """
    
    :param sc: 
    :param session: 
    :return: 
    """
    embedDict = initEmbedDict(session)
    prepareCoordinates(embedDict)
    TForm = getFieldTFormPar(sc, embedDict, zExpand=False, xDim=(2, 70))
    embedDict['TForm1'] = TForm
    session.embedDict = embedDict
    if return_t:
        return TForm
    else:
        return session


def getTFromExp(sc, session):
    """
    
    :param sc: 
    :param session: 
    :return: 
    """
    embedDict = initEmbedDict(session)
    prepareCoordinatesZExp(embedDict, session.regDict['grpZPos'], True)
    return getFieldTFormPar(sc, embedDict, zExpand=True, xDim=(2, 70))


def showEmbedding(img, TForm, sc=None, figure_size=(20, 20)):
    """

    :param img: image to show
    :param TForm: transform to show
    :param sc: Spark Context
    :param figure_size: figure size to plot
    :return: plots the image and the number of the field in the upper left corner
    """
    if len(img.shape) == 3:
        if sc is not None:
            img = getExampleVol(sc, img, TForm, project=True)
        else:
            logger.error('3d img given with no sc')
            return
    plt.figure(figsize=figure_size)
    plt.imshow(img)
    locations = np.min(TForm['brainFields'], axis=(0, 1))[:2, :]
    for i, (x, y) in enumerate(locations.T):
        plt.text(y, x, i, fontsize=22, color='red')
