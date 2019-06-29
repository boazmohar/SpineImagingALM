# Encoding: utf-8
""" HeaderReader class for getting tif header information saved by ScanImage

"""

from __future__ import print_function
import re
import copy
import numpy as np


class HeaderReader(object):
    """Loader object used to read tiff headers
    """

    def __init__(self, sparkContext):
        """ Initialize a new HeaderReader object.

        :param sparkContext: The pyspark SparkContext object used by the current Thunder environment.
        """
        self.sc = sparkContext

    def read(self, dataPath, planesPerVol):
        """ parallel reading of tif files to ger header information

        :param dataPath: folder to tiff files
        :param planesPerVol: number of planes per volume
        :return: stamps a n by 4 array with: [0] - frame number; [1] - frame time stamp
                 [2] - Acq Trigger Timestamps; [3] Next File Marker Timestamps
        """
        def _getStamps(myData):

            from skimage.external import tifffile
            filename = myData[0]
            nPlanes = myData[1]

            tfh = tifffile.TiffFile(filename)
            nPages = len(tfh.pages)
            if nPages % nPlanes != 0:
                print('error: incomplete volume')
                # return None
            customHeaders = np.zeros((nPlanes, 4))
            out = []
            for page in range(0, nPages):
                k = (page + 1) % nPlanes
                if k == 0:
                    k = nPlanes

                desStr = tfh.pages[page].tags['image_description'].value.decode('utf-8')
                desStr = " ".join(desStr.split()).replace('\\', ' ')
                desStr = desStr.replace(')', '')
                desStr = desStr.replace('(', '')
                try:
                    customHeaders[k - 1, 0] = re.search(
                        '(?<=Frame\sNumber\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                    customHeaders[k - 1, 1] = re.search(
                        '(?<=Frame\sTimestamps\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                    customHeaders[k - 1, 2] = re.search(
                        '(?<=Acq\sTrigger\sTimestamps\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                    customHeaders[k - 1, 3] = re.search(
                        '(?<=Next\sFile\sMarker\sTimestamps\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                except AttributeError:
                    #SI2016
                    customHeaders[k - 1, 0] = re.search(
                        '(?<=frameNumbers\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                    customHeaders[k - 1, 1] = re.search(
                        '(?<=frameTimestamps_sec\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                    customHeaders[k - 1, 2] = re.search(
                        '(?<=acqTriggerTimestamps_sec\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                    customHeaders[k - 1, 3] = re.search(
                        '(?<=nextFileMarkerTimestamps_sec\s=\s)[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?',
                        desStr).group(0)
                if k == nPlanes:
                    out.append(copy.deepcopy(customHeaders))
            return out

        import glob
        files = glob.glob(dataPath + '*.tif')
        data = [None] * len(files)
        for v in range(0, len(files)):
            data[v] = (files[v], planesPerVol)
        stamps = self.sc.parallelize(data).flatMap(_getStamps).collect()
        stamps = sorted(stamps, key=lambda stamp: stamp[0, 0])
        return stamps
