# Encoding: utf-8
""" General utilities disk IO.

"""
from __future__ import print_function, division

import copy
import logging
import os
import shelve
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from pySparkUtils.utils import balanced_repartition, save_rdd_as_pickle, load_rdd_from_pickle
from skimage.external import tifffile

from prep.Embedding import getExampleVol

# Setup logging
logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def cropEmbed(sc, data, path, dims, TForm, name, writeFlag=True, returnFlag=False):
    """ crop a 3d or 4d array, embed using TForm (projected) and write to path given or return

    :param sc: Spark context
    :param data: 3d or 4d array
    :param path: session view folder usually
    :param dims: dimensions of the original volume before expanding
    :param TForm: TForm dict
    :param name: name of file to write
    :param writeFlag: if True will write file to disk
    :param returnFlag if True will return the resulting embedding
    :return: depending on return flag
    """
    sz = data.shape
    if len(sz) == 4:
        sz = sz[1:]
    elif len(sz) != 3:
        raise ValueError('Data must be 3 or 4 dimentions got %d' % len(sz))
    temp = copy.deepcopy(data)
    x_diff = float(sz[0] - dims[0]) / 2.0
    x_start = int(np.floor(x_diff))
    x_stop = np.ceil(x_diff)
    x_stop = int(sz[0] - x_stop)

    y_diff = float(sz[1] - dims[1]) / 2.0
    y_start = int(np.floor(y_diff))
    y_stop = np.ceil(y_diff)
    y_stop = int(sz[1] - y_stop)

    if len(data.shape) == 4:
        temp = temp[:, x_start:x_stop, y_start:y_stop, :]
    else:
        temp = temp[x_start:x_stop, y_start:y_stop, :]
    tempEmbed = getExampleVol(sc, temp, TForm, project=True)
    if writeFlag:
        if len(data.shape) == 4:
            writeTiff(path, tempEmbed.transpose((1, 2, 0)), name)
        else:
            writeTiff(path, np.expand_dims(tempEmbed, 2), name)
    if returnFlag:
        return tempEmbed


def writeTiff(path, image, name, dtype='float32'):
    """ saves a 2d or 3d numpy array to multi-page tiff file as (x,y,z)

    :param path: base path
    :param image: 2d / 3d numpy array
    :param name: image name (.tif will be added)
    :param dtype: dtype of image (float32)
    """
    filename = path + name + ".tif"
    try:
        if image.size * image.dtype.itemsize > 2000 * 2 ** 20:
            writer1 = tifffile.TiffWriter(filename, bigtiff=True)
        else:
            writer1 = tifffile.TiffWriter(filename)
        if len(image.shape) > 2:
            writer1.save(image.transpose((2, 0, 1)).astype(dtype=dtype), photometric="minisblack")
        else:
            writer1.save(image.astype(dtype=dtype), photometric="minisblack")
        writer1.close()
        del writer1
        os.chmod(filename, 0o777)
    except (IOError, OSError):
        myDir, myName = os.path.split(filename)
        newFilename = tempfile.mktemp(suffix='_' + myName, prefix='tmp', dir=myDir)
        if image.size * image.dtype.itemsize > 2000 * 2 ** 20:
            writer1 = tifffile.TiffWriter(filename, bigtiff=True)
        else:
            writer1 = tifffile.TiffWriter(filename)
        logger.error('error saving tiff:' + str(filename) + ' changed name to temp: ', str(newFilename))
        if len(image.shape) > 2:
            writer1.save(image.transpose((2, 0, 1)).astype(dtype=dtype), photometric="minisblack")
        else:
            writer1.save(image.astype(dtype=dtype), photometric="minisblack")
        writer1.close()
        del writer1
        os.chmod(newFilename, 0o777)


def writeTiff4D(path, image, name, dtype='float32'):
    """ saves a 4d numpy array to multi-page tiff as (x,y,z*t)

    :param path: base path
    :param image: 4d numpy array assumes 2nd and 3rd dimensions are x,y
    :param name: image name (.tif will be added)
    :param dtype: dtype of image (float32)
    """
    image = image.transpose((3, 0, 1, 2))
    si = image.shape
    filename = path + name + ".tif"
    if image.size * image.dtype.itemsize > 2000 * 2 ** 20:
        writer1 = tifffile.TiffWriter(filename, bigtiff=True)
    else:
        writer1 = tifffile.TiffWriter(filename)
    out = image.reshape(si[0] * si[1], si[2], si[3], order="F")
    writer1.save(out.astype(dtype=dtype), photometric="minisblack")
    writer1.close()
    del writer1
    os.chmod(path + name + ".tif", 0o777)


def loadmat(filename):
    """    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    :param filename: filename to be loaded
    :return loaded .mat file as dictionary
    """

    def _check_keys(dict_mat):
        """
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in dict_mat:
            if isinstance(dict_mat[key], sio.matlab.mio5_params.mat_struct):
                dict_mat[key] = _todict(dict_mat[key])
        return dict_mat

    def _todict(matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        dict_mat = {}
        for string in matobj._fieldnames:
            elem = matobj.__dict__[string]
            if isinstance(elem, sio.matlab.mio5_params.mat_struct):
                dict_mat[string] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                dict_mat[string] = _tolist(elem)
            else:
                dict_mat[string] = elem
        return dict_mat

    def _tolist(array):
        """
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        elem_list = []
        for sub_elem in array:
            if isinstance(sub_elem, sio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list

    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def loadData(sc, session, cutoff=None, saveBinary=True, start=None, stop=None, xStart=2, xStop=70, overwrite=False,
             timepoints=None, repartition=True, cutoff_fallback=40, zoom=None, return_clean_path=False,
             binary_path=None, **kwargs):
    """ Loads raw data, and return clean after cropping using a cutoff to get rig of electrical noise and correcting for
    bi-directional phase offset

    :param sc: Spark Context
    :param session: SpineSession object
    :param cutoff: cutoff value. if None will try to infer from flyback frames
    :param saveBinary: if True will save to nrs as binary with one file per partition
    :param start: start file to load
    :param stop: last file to load
    :param xStart: first rows to crop
    :param xStop: last rows to crop
    :param overwrite: if true will overwrite the binary saved
    :param timepoints: how many timepoint to use to estimate the bi-directional phase offset
    :param repartition: if True will repartition the data into sc.defaultParallelism * 2 for better performance
    :param cutoff_fallback: if cutoff is None and there are no flyback frames will use this value
    :param zoom: zoom
    :param binary_path: base path to save binary file
    :param kwargs: to be passed to initBase
    :return: raw data Images object, clean Images object
    """
    session.initBase(sc, **kwargs)
    data = session.loadRawData(sc=sc, start=start, stop=stop)
    if repartition:
        data = balanced_repartition(data, sc.defaultParallelism * 2)
    if cutoff is None:
        if not hasattr(session, 'fieldMaskBool') or len(np.where(session.fieldMaskBool == 0)[0]) == 0:
            cutoff = cutoff_fallback
            logger.info('Setting cutoff at %d, No flyback frames' % cutoff_fallback)
        else:
            data.cache()
            data.count()
            flyback = np.where(session.fieldMaskBool == 0)[0]
            flybackData = data[:, :, :, flyback]
            flybackDataMean = flybackData.map(np.mean).toarray()
            flybackDataSTD = flybackData.map(np.std).toarray()
            m = np.median(flybackDataMean).astype(float)
            s = np.median(flybackDataSTD).astype(float)
            cutoff = m + 2.5 * s
            logger.info('Mean: %.2f, STD: %.2f, Cutoff: %.2f' % (m, s, cutoff))
    session.display()
    plt.show()
    if timepoints is None:
        timepoints = data.shape[0]
    session.addPipeline('Clean')
    session.addStep('Clean', 'crop1', 'cropDataStep', planes=session.fieldMask)
    session.addStep('Clean', 'clean', 'cleanDataStep', cutoff=cutoff)
    if zoom is not None:
        session.addStep('Clean', 'zoom', 'zoomDataStep', zoom=zoom)
    session.addStep('Clean', 'scanPhase', 'scanPhaseStep', timepoints=timepoints)
    session.addStep('Clean', 'crop2', 'cropDataStep', xStart=xStart, xStop=xStop)
    clean = session.doPipeline(pipeline='Clean', data=data)
    clean.cache()
    clean.count()
    data.uncache()
    if binary_path is None:
        name = '/nrs/svoboda/moharb/New/'
    else:
        name = binary_path
    name = name + session.animalID + '_' + session.date + session.run + 'CleanBinaryPickle'
    if saveBinary:
        save_rdd_as_pickle(clean, name, overwrite=overwrite)
    if return_clean_path:
        return data, clean, name
    else:
        return data, clean


def reloadClean(sc, session, name=None, returnRegDict=False, returnTForm=False, repartition=True, full_path=None):
    """ reloads session and clean data from nrs

    :param sc: Spark context
    :param session: SpineSession object
    :param name: name of session to load
    :param returnRegDict: if to return regDict
    :param returnTForm:  if to return TForm
    :param repartition: if to repartition Clean RDD. True = sc.defaultParallelism * 2 or any other int > 0
    :param full_path: if the full path of the data is known (other then cluster default)
    :return: either clean (if all is set to None / False) or: clean, session, regDict, TForm1
    """
    if name is not None:
        session = session.load(name)
    if full_path is None:
        full_path = '/nrs/svoboda/moharb/New/' + session.animalID + '_' + session.date + session.run + 'CleanBinaryPickle'
    clean = load_rdd_from_pickle(sc, full_path)
    # if bool true
    if repartition > 0:
        # if number specified
        if repartition > 1:
            nPartitions = repartition
        # if number not specified
        else:
            nPartitions = sc.defaultParallelism * 2
        clean = clean.repartition(nPartitions)
        logger.info('Repartitioned Clean to %d partitions' % nPartitions)
    clean.cache()
    clean.count()
    if returnRegDict:
        regDict = session.regDict
        regDict['data'] = clean
    else:
        regDict = None
    if returnTForm:
        if hasattr(session, 'embedDict') and 'TForm1' in session.embedDict:
            TForm1 = session.embedDict['TForm1']
        else:
            logger.error('TForm Not Found')
            TForm1 = None
    else:
        TForm1 = None
    if name is None and regDict is None and TForm1 is None:
        return clean
    else:
        return clean, session, regDict, TForm1


def saveRegData(session, regData, overwrite=False, base='/nrs/svoboda/moharb/New/', return_name=False):
    """ saves registered data Images object

    :param session: SpineSession object
    :param regData: thunder Images object
    :param overwrite: if True will overwrite files
    :param base: base path to save to
    """
    name = base + session.animalID + '_' + session.date + session.run + 'RegBinaryPickle'
    save_rdd_as_pickle(regData, name, overwrite=overwrite)
    if return_name:
        return name


def loadRegData(sc, session, name='step5', base='/nrs/svoboda/moharb/New/', check=True):
    """ loads regData and session if name is not None

    :param sc: spark context
    :param session: session object
    :param name: session name. If None will not load a session
    :param base: base path of sessions
    :return: session, regData
    """
    if name is not None:
        session = session.load(name)
    path = base + session.animalID + '_' + session.date + session.run + 'RegBinaryPickle'
    regData = load_rdd_from_pickle(sc, path)
    regData.cache()
    count = regData.count()
    if check and count != session.regDict['globalTC'].shape[0]:
        logger.error('count: %d different from globalTC length %d' % (count, session.regDict['globalTC'].shape[0]))
    logger.info('Loaded regData from: %s' % path)
    return session, regData


def save_shelve(variables, base='/groups/svoboda/svobodalab/users/Aaron/Database/Vision/', filename='shelve',
                prefix=None, flag='c', verbose=True):
    """

    :param base: folder name
    :param filename: file name
    :param prefix: only save variables starting with this prefix (could be a tuple)
    :param flag: 'n' Always create a new, empty database, open for reading and writing
                 'w' Open existing database for reading and writing
                 'c' Open database for reading and writing, creating it if it doesn’t exis
    :return:
    """
    filename_full = os.path.join(base, filename)
    my_shelf = shelve.open(filename_full, flag)
    for key in variables:
        if prefix is None or key.startswith(prefix):
            try:
                my_shelf[key] = globals()[key]
                if verbose:
                    logger.info('Saved: %s' % key)
            except Exception:
                logger.error('ERROR shelving: {0}'.format(key))


def load_shelve(base='/groups/svoboda/svobodalab/users/Aaron/Database/Vision/', filename='shelve', prefix=None):
    """

    :param base: folder name
    :param filename: file name
    :param prefix: only load variables starting with this prefix (could be a tuple)
    :return:
    """
    filename_full = os.path.join(base, filename)
    my_shelf = shelve.open(filename_full)
    for key in my_shelf:
        if prefix is None or key.startswith(prefix):
            globals()[key] = my_shelf[key]
    my_shelf.close()
