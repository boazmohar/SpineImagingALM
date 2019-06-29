# Encoding: utf-8
""" Session class
 Session is used to handle imaging data from ScanImage using Thunder and Spark
 It has load and save capabilities

"""
from __future__ import print_function

import copy
import glob
import os
import sys
import logging

import matplotlib
import pyspark
import thunder as td

from .IO import writeTiff, writeTiff4D
from .Pipeline import Pipeline
from .Utils import getTiffInfo, searchTiffInfo

try:
    import cPickle as pickle
    input = raw_input
    PY3 = False
except ImportError:
    import pickle
    PY3 = True

logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)

class Session(object):
    """Session object used to analyze imaging data from ScanImage
    """

    def __init__(self, basePath='/groups/svoboda/svobodalab/users/Aaron', animalID=None, date=None, run=None,
                 nPlanes=None, verbose=True):
        """ Initialize a new Session object. Assumes a format of BASE\YYMMDD_ID\RUN
        :param basePath: base data folder
        :param animalID: such as WR34
        :param date: as YYMMDD format
        :param run: name of run
        :param nPlanes: number of planes
        """
        self.animalID = animalID
        self.date = date
        self.run = run
        self.basePath = basePath
        self.nPlanes = nPlanes
        if basePath and animalID and date and run:
            self.path = os.path.join(basePath, date + '_' + animalID, run)
            self.viewPath = os.path.join(self.path + 'View', '')
            try:
                if not os.path.isdir(self.viewPath):
                    os.mkdir(self.viewPath)
            except OSError:
                logger.error('Could not make directory: %s' % self.viewPath)
                self.viewPath = None
        elif verbose:
            logger.info('Either basePath, animalID, date or run were missing')
        self.pipelines = dict()

    def __repr__(self):
        if self.animalID and self.date and self.run:
            return 'Session: animal: %s, date: %s, run: %s' % (self.animalID, self.date, self.run)
        else:
            return 'Session object'

    def display(self):
        """ method to print session info  """
        print('Session object with properties: ')
        prop_to_print = ['path', 'viewPath', 'volRate']
        for prop in prop_to_print:
            if hasattr(self, prop):
                print(str(prop) + ' = ' + str(getattr(self, prop)))

    def initBase(self, sc=None, nPartitions=None, flyVelocity=266):
        """ initialize properties of the session

        :param sc: SparkContext
        :param nPartitions: number of partitions to load the data, if None will take  sc.defaultParallelism
        :param flyVelocity: X, Y galvos speed in um/ms
        """

        self.flyVelocity = flyVelocity

        # getting tiff header information
        try:
            self.tiffHeader = getTiffInfo(self.path)
            self.volRate = float(searchTiffInfo(self.tiffHeader, 'scanVolumeRate'))
        except (AttributeError, IndexError):
            logger.error('No tiff header found')

        if nPartitions:
            self.nPartitions = nPartitions
        elif sc is not None:
            self.nPartitions = sc.defaultParallelism
        else:
            logger.error('Can not determine number of partitions, setting to 1')
            self.nPartitions = 1

    def loadRawData(self, sc, start=None, stop=None):
        """ loads data from disk using thunder

        :param sc: spark context to load images with
        :param start: file number to start loading
        :param stop: file number to end loading
        :return: thunder images object
        """
        # loading images
        if start is None:
            start = 0
        self.start = start
        if stop is None:
            files = glob.glob(os.path.join(self.path, '') + '*.tif')
            stop = len(files)
        self.stop = stop
        data = td.images.fromtif(self.path, nplanes=self.nPlanes, engine=sc, npartitions=self.nPartitions,
                                 start=start, stop=stop, discard_extra=True)
        if self.nPlanes > 1:
            data = data.map(lambda x: x.transpose(1, 2, 0))
        data.cache()
        logger.info('shape: ' + str(data.shape))
        self.data_shape = data.shape
        return data

    def load(self, name):
        """ reassign your spineSession object with one previously saved to disk
        use as: session = session.load('myName')
        :param name:  name of pickle file to be loaded
        :return: a new spine session object
        """
        with open(self.path + name + '.p', 'rb') as inputData:
            if PY3:
                a = pickle.load(inputData, encoding='latin1')
            else:
                a = pickle.load(inputData)
            logger.info('Session loaded from ' + self.path + name)
            return a

    def save(self, name, confirm=True):
        """ saves a Session object to disk, will remove any references to RDDs
        saving RDDs should be done separately
        :param name: name to give the Pickle file
        :param confirm: whether to ask if to overwrite the current file
        """

        def searchDict(myDict):
            rddBackup = {}
            for key in list(myDict.keys()):
                if isinstance(myDict[key], (td.base.Base, pyspark.RDD, pyspark.Broadcast)):
                    rddBackup[key] = myDict[key]
                    del myDict[key]
                elif isinstance(myDict[key], matplotlib.figure.Figure):
                    del myDict[key]
            myDict2 = copy.deepcopy(myDict)
            for key in list(rddBackup.keys()):
                myDict[key] = rddBackup[key]
            return myDict2

        Vars = vars(self).keys()
        newSession = self.__class__(verbose=False)
        for var in Vars:
            current = getattr(self, var)
            if isinstance(current, dict) and not var == 'Sp':
                setattr(newSession, var, searchDict(getattr(self, var)))
            # if it is not an RDD or matplotlib object
            elif not isinstance(current, (td.base.Base, pyspark.RDD, matplotlib.figure.Figure,
                                          logging.Logger)):
                setattr(newSession, var, getattr(self, var))
        filename = self.path + name + '.p'
        if os.path.isfile(filename) and confirm:
            answer = input('File \'' + filename + '\' exists\noverwrite? (y-Yes, n-No): ')
        else:
            answer = 'y'
        if answer == 'y':
            with open(self.path + name + '.p', 'wb') as output:
                pickle.dump(newSession, output, 2)
            logger.info('Session saved to  \'' + filename + '\'')

    def writeTiff(self, image, name, path=None, type2='float32'):
        """ wrapper function with default path to session.viewPath """
        if not path:
            path = self.viewPath
        writeTiff(path, image, name, type2)

    def writeTiff4D(self, image, name, path=None, type2='float32'):
        """ wrapper function with default path to session.viewPath """
        if not path:
            path = self.viewPath
        writeTiff4D(path, image, name, type2)

    def findProp(self, name, default=None):
        """ finds a property, also looks inside dictionaries

        :param name: name of property to look for
        :param default: what to return if didn't find it
        :return: property value or default if not found
        """
        Vars = vars(self).keys()
        for var in Vars:
            current = getattr(self, var)
            if isinstance(current, dict) and not var == 'Sp':
                myDict = getattr(self, var)
                if name in myDict.keys():
                    return myDict[name]
            elif var == name:
                return getattr(self, var)
        return default

    def addPipeline(self, pipeline, source='raw'):
        self.pipelines[pipeline] = Pipeline(source)

    def delPipeline(self, pipeline):
        if pipeline in self.pipelines:
            del self.pipelines[pipeline]

    def addStep(self, pipeline, name, func, **kwargs):
        self.pipelines[pipeline].addStep(name, func, **kwargs)

    def delStep(self, pipeline, name):
        if pipeline in self.pipelines:
            self.pipelines[pipeline].remStep(name)
        else:
            logger.error('No pipeline named:' + str(pipeline))

    def doPipeline(self, pipeline, data):
        if pipeline in self.pipelines:
            data = self.pipelines[pipeline].doAllSteps(data)
            return data
        else:
            logger.error('No pipeline named:' + str(pipeline))
