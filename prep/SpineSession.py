# Encoding: utf-8
""" SpineSession class
 Spine session is used to load metadata information created using matlab SpineImaging class
 It has load and save capabilities
"""
import copy
import io
import json
import os
import sys
import pickle
import logging
from builtins import input

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyspark
import thunder as td
import requests

from .IO import loadmat
from .Session import Session
from .Utils import searchTiffInfo

logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)


class SpineSession(Session):
    """Session object used to analyze spine imaging data
    """

    def __init__(self, basePath='/groups/svoboda/svobodalab/users/Aaron', animalID=None, date=None, run=None,
                 verbose=True, nPlanes=0):
        """ Initialize a new SpineSession object.
        Assumes a format of BASE\YYMMDD_animalID\RUN
        If no input is given, returns an empty object
        :param basePath: where the base folder is in
        :param animalID: such as WR34
        :param date: as YYMMDD format
        :param run: name of run
        :param verbose: Use logger
        :param nPlanes number of planes
        """
        super(SpineSession, self).__init__(basePath, animalID, date, run, nPlanes=nPlanes, verbose=verbose)
        if animalID and date and run:
            logger.info('SpineSession initialized, path = ' + self.path)

    def __repr__(self):
        if self.animalID and self.date and self.run:
            return 'SpineSession: animal: %s, date: %s, run: %s' % (self.animalID, self.date, self.run)
        else:
            return 'SpineSession object'

    def display(self):
        """ method to print session info  """
        super(SpineSession, self).display()
        if hasattr(self, 'ZError'):
            Last = self.TrueFieldsCenterSamp[-1]
            Errors = self.ZError[self.TrueFieldsCenterSamp]

            fig = plt.figure()
            plt.subplot(2, 1, 1)
            plt.plot(self.ZError[0:Last + 10])
            plt.plot(self.TrueFieldsCenterSamp, Errors, 'o')
            plt.legend(['Raw error', 'Scan fields'])
            plt.ylabel('Error (um)')
            plt.xlabel('Lines')
            plt.title('Z errors')
            self.ZErrorFig = fig

    def initBase(self, sc=None, xScale=1.13, yScale=1.13, nPartitions=None, xScaleAnat=1.0, yScaleAnat=1.0,
                 baseFlyLines=4, pixSizeZ=2500, flyVelocity=266):
        """  initialize properties of the session

        :param sc: SparkContext
        :param xScale: expansion factor for pixel in X
        :param yScale: expansion factor for pixel in X
        :param nPartitions: number of partitions to load the data, if None will take  sc.defaultParallelism
        :param flyVelocity: X, Y galvos speed in um/ms
        :param xScaleAnat: scaling factor for stack in X
        :param yScaleAnat: scaling factor for stack in Y
        :param baseFlyLines: number of fly lines to remove in any case
        :param pixSizeZ: PSF FWHM in Z in nm
        """

        super(SpineSession, self).initBase(sc, nPartitions, flyVelocity)
        self.xScale = xScale
        self.yScale = yScale
        self.xScaleAnat = xScaleAnat
        self.yScaleAnat = yScaleAnat
        self.pixSizeZ = pixSizeZ
        self.baseFlyLines = baseFlyLines
        self.getSpMat()
        if hasattr(self, 'xSize'):
            self.pixSizeXY = 1.0 / (self.xSize / self.xSizeOrig) * xScale * 1000  # in nm!!!
            self.setFlyLines(self.baseFlyLines)
        else:
            self.pixSizeXY = None
        self.getMeta(from_web=True)

    def getMeta(self, filename='/Database/Sessions.csv', from_web=False,
                key='1ENbQe4QtUZ6as9y8Kcr_sHl3YpwQoAs_EScIIURefSQ', gid='78168970'):
        """

        :param from_web: if true will go to public google sheet, if false will go to session.basePath + filename
        :param filename: local file to look for
        :param key: google sheet key
        :param gid: google sheet gid
        :return: metadata from the session
        """
        if from_web:
            response = requests.get('https://docs.google.com/spreadsheet/ccc?key=' + key + '&output=csv&gid=' + gid)
            assert response.status_code == 200, 'Wrong status code'
            f = io.StringIO(response.content.decode('utf-8'))
            sessionsDF = pd.read_csv(f)
        else:
            # read the file
            databasePath = self.basePath + filename
            sessionsDF = pd.read_csv(databasePath)
        # get the animal ID
        index1 = sessionsDF.WR == self.animalID
        # get the date
        startTime = sessionsDF.StartTime
        startTime = startTime.apply(lambda x: pd.to_datetime(x).date())
        index2 = startTime == pd.to_datetime(self.date, yearfirst=True).date()
        # get the run
        run = int(self.run[3])
        index3 = sessionsDF.Run == run
        index4 = sessionsDF.Type == 'MROI'
        # find the session
        self.meta = sessionsDF[index1 & index2 & index3 & index4].to_dict('list')
        json.dump(self.meta, open(self.path + 'meta.txt', 'w'))

    def getSpMat(self):
        """ loads the information from the Sp.mat file of this session

        :return: updated SpineSession object
        """
        full_path = os.path.join(self.path, 'Sp.mat')
        if os.path.isfile(full_path):
            Sp = loadmat(full_path)
            Sp = Sp['Sp_s']
            self.Sp = Sp
            self.nPlanes = len(Sp['All']['x'])
            fieldMask = Sp['TrueFieldsMask']
            Options = Sp['OptionsStruct']
            self.lineRate = Options['lineRate']
            self.optFlyLines = Options['FlyLines']
            ySize = Options['ImagingLines']
            self.ySize = ySize
            xSize = Options['ImagingPixels']
            self.xSize = xSize
            self.fieldMask = fieldMask
            self.fieldMaskBool = np.bool_(fieldMask)
            self.rawX = np.array(Sp['All']['x'])
            self.rawY = np.array(Sp['All']['y'])
            self.rawZ = np.array(Sp['All']['z'])
            self.x = self.rawX[self.fieldMaskBool]
            self.y = self.rawY[self.fieldMaskBool]
            self.z = self.rawZ[self.fieldMaskBool]
            self.xSizeOrig = Options['xSize']
            self.ySizeOrig = Options['ySize']
            self.TrueFieldsCenterSamp = Sp['TrueFieldsCenterSamp']
            self.ZWave = Sp['Zwave']
            self.UM_1X = Options['UM_1X']
            # Z error loading
            if 'ZError' in Sp:
                self.ZError = np.array(Sp['ZError']).flatten()
                self.ZActual = self.ZWave + self.ZError
                logger.info('getSpMat:: Loaded ZError')
            else:
                logger.error('getSpMat:: No ZError found')
        else:
            logger.error('No Sp.mat file found in : %s' % full_path)

    def setFlyLines(self, baseLines):
        """ set fly lines to minimum baseLines minus number of flyto lines from SI

        :param baseLines: minimum lines to remove not excluding SI flyto
        :return: updates self.flyLines
        """

        # look for hResScan (SI2015) or hScan2D (SI2016 which moved to BigTiff)
        try:
            val = searchTiffInfo(self.tiffHeader, 'hResScan.flytoTimePerScanfield')
        except AttributeError:
            logger.info('SI2016 Big Tiff')
            # SI2016
            val = searchTiffInfo(self.tiffHeader, 'hScan2D.flytoTimePerScanfield')
        self.baseFly = baseLines - round(float(val) * self.lineRate)
        if self.fieldMaskBool[-1]:
            # if there is no flyback field!
            X = copy.deepcopy(np.append(self.rawX, self.rawX[0]))
            Y = copy.deepcopy(np.append(self.rawY, self.rawY[0]))
            flyDist = np.max(np.concatenate((np.absolute(np.diff(X).reshape(-1, 1)),
                                             np.absolute(np.diff(Y).reshape(-1, 1))), axis=1), axis=1)
            flyDist = flyDist[self.fieldMaskBool]
        else:
            X = copy.deepcopy(self.rawX)
            Y = copy.deepcopy(self.rawY)
            flyDist = np.max(np.concatenate((np.absolute(np.diff(X).reshape(-1, 1)),
                                             np.absolute(np.diff(Y).reshape(-1, 1))), axis=1), axis=1)
            flyDist = flyDist[self.fieldMaskBool[:-1]]
        flyLines = np.round(flyDist / self.flyVelocity * self.lineRate / float(1000)).astype('int64')
        flyLines[flyLines < self.baseFly] = self.baseFly
        self.flyLines = flyLines

    def showLocations(self):
        """

        :return:
        """
        zStep = self.Sp['OptionsStruct']['Zstep']
        refPx = self.Sp['OptionsStruct']['RefPixels']
        UM_1X = self.Sp['OptionsStruct']['UM_1X']
        x = self.Sp['All']['x']
        y = self.Sp['All']['y']
        z = self.Sp['All']['z']
        xyStep = float(UM_1X) / float(refPx)
        for i in range(len(np.where(self.fieldMask)[0])):
            index = np.where(self.fieldMask)[0][i]
            logger.info('Cell:%d, X:%.2f, Y:%.2f, Z:%d' % (i,
                                                           x[index] / xyStep,
                                                           y[index] / xyStep,
                                                           z[index] / zStep))
