# Encoding: utf-8
""" Behavior module for loading data from RT linux behavior controller.

"""

from __future__ import print_function

import copy
import matplotlib.pyplot as plt
import logging
import numpy as np
import os
import sys

from prep.Embedding import getExampleVol
from prep.HeaderReader import HeaderReader
from prep.IO import loadmat

# Setup logging
logger = logging.getLogger(__name__)
logger.handlers = []
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter("%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"))
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def loadBehavior(sc, session, filename='Behavior'):
    """ Loads behavior data from RT linux mat file and the timestamps from
    the imaging files

    :param sc: SparkContext
    :param session: current session object
    :param filename: file name
    :return behavDict: Behavior dictionary
    """
    behavDict = dict()
    behavDict['filename'] = filename
    behavDict['path'] = session.path
    behavDict['nPlanes'] = session.nPlanes
    if hasattr(session, 'start') and session.start is not None:
        behavDict['start'] = session.start
    else:
        behavDict['start'] = 1

    # from RT linux
    fullPath = os.path.join(behavDict['path'], behavDict['filename'])
    BehavData = loadmat(fullPath)
    behavDict['rawData'] = BehavData
    Saved = BehavData['saved']
    behavDict['Delay'] = Saved['TimesSection_Delay_period']
    if hasattr(session, 'stop') and session.stop is not None:
        behavDict['stop'] = session.stop
    else:
        behavDict['stop'] = len(behavDict['rawData']['saved_history']['TimesSection_title'])
    try:
        behavDict['CorrectR'] = Saved['yes_no_multi_pole_akobj_correct_R_history']
        behavDict['CorrectL'] = Saved['yes_no_multi_pole_akobj_correct_L_history']
        behavDict['IncorrectR'] = Saved['yes_no_multi_pole_akobj_incorrect_R_history']
        behavDict['IncorrectL'] = Saved['yes_no_multi_pole_akobj_incorrect_L_history']
        behavDict['EarlyLick'] = Saved['yes_no_multi_pole_akobj_early_history']
    except KeyError:
        behavDict['CorrectR'] = Saved['yes_no_multi_pole_switchobj_correct_R_history']
        behavDict['CorrectL'] = Saved['yes_no_multi_pole_switchobj_correct_L_history']
        behavDict['IncorrectR'] = Saved['yes_no_multi_pole_switchobj_incorrect_R_history']
        behavDict['IncorrectL'] = Saved['yes_no_multi_pole_switchobj_incorrect_L_history']
        behavDict['EarlyLick'] = Saved['yes_no_multi_pole_switchobj_early_history']

    # verify consistency with image timestamps
    History = BehavData['saved_history']
    Events = History['RewardsSection_LastTrialEvents']
    Events = Events[1:]  # disregard first trial!!!
    behaviorTime = np.zeros(len(Events))
    behavDict['behaviorTime'] = behaviorTime

    # get the transition to state 40 from RT linux as the trigger time
    i = 0
    for event in Events:
        FirstArray = np.array([np.array(xi) for xi in event])
        Index40 = np.where(FirstArray[:, 0] == 40)[0][1]
        behaviorTime[i] = FirstArray[Index40, 2]
        i += 1
    behaviorTime = np.diff(behaviorTime)

    # get the timestamps from tiff header
    hrObj = HeaderReader(sc)
    stamps = hrObj.read(os.path.join(behavDict['path'], ''), behavDict['nPlanes'])
    behavDict['stamps'] = stamps
    behavDict['nRecords'] = len(stamps)

    # check missing frames
    stamps_frames = np.asarray(stamps)
    sz = stamps_frames.shape
    stamps_frames = np.reshape(stamps_frames, (sz[0] * sz[1], sz[2]))
    behavDict['stamps_frames'] = stamps_frames
    if np.sum(np.diff(stamps_frames[:, 0]) - 1) > 0:
        logger.info('loadBehavior:: missing frames')
    else:
        logger.info('loadBehavior:: no missing frames')

    # get the time where the timestamp of the trial changed
    ImageTimes = np.diff(stamps_frames[np.where(np.diff(stamps_frames[:, 3]) > 0), 3])
    ImageTimes = np.squeeze(ImageTimes.T)
    behavDict['imageTime'] = ImageTimes

    # plot times
    plt.plot(behaviorTime, label='Behavior')
    plt.plot(ImageTimes, '*', label='Imaging')
    plt.legend()
    plt.xlabel('Trial #')
    plt.ylabel('Trial duration (s)')

    # Num of Correct trials:
    CR = np.array(behavDict['CorrectR'], dtype=np.bool)
    CL = np.array(behavDict['CorrectL'], dtype=np.bool)
    EL = np.array(behavDict['EarlyLick'], dtype=np.bool)
    numCR = sum(np.logical_and(CR, np.logical_not(EL)))
    numCL = sum(np.logical_and(CL, np.logical_not(EL)))
    total = sum(np.logical_not(EL))
    behavDict['numCR'] = numCR
    behavDict['numCL'] = numCL
    behavDict['numER'] = sum(np.logical_and(behavDict['IncorrectR'], np.logical_not(behavDict['EarlyLick'])))
    behavDict['numEL'] = sum(np.logical_and(behavDict['IncorrectL'], np.logical_not(behavDict['EarlyLick'])))
    behavDict['total'] = total
    CR = Saved['AnalysisSection_PercentCorrect_R']
    CL = Saved['AnalysisSection_PercentCorrect_L']
    behavDict['CR'] = CR
    behavDict['CL'] = CL
    plt.title('CR: %d %f%%, CL: %d %f%%, ER:, %d, EL:, %d' % (numCR, CR, numCL, CL, behavDict['numER'],
                                                              behavDict['numEL']))

    # get two level index
    getIndices(behavDict)
    return behavDict


def getIndices(behavDict):
    """ prepares the index for aggregation of the data according to behavior

    :param behavDict: Behavior dictionary after loadBehavior runs
    """

    nRecords = behavDict['nRecords']
    StartTrial = np.asarray(np.where(np.diff(behavDict['stamps_frames'][:, 3]) > 0))
    StartTrial = np.concatenate((np.asarray([[0]]), StartTrial), axis=1)
    StartTrial = StartTrial[0]
    StartTrial = StartTrial / behavDict['nPlanes']
    Index = np.zeros(nRecords, dtype=np.int)
    for i in StartTrial:
        Index[int(i)::] = Index[int(i)::] + 1
    if behavDict['stop'] is None:
        behavDict['stop'] = np.max(Index)
    TimeIndex = np.zeros(nRecords, dtype=np.int)  # frame number from start of current trial
    TypeIndex = np.zeros(nRecords, dtype=np.int)  # trial type: 1 CR, 2 CL, 3 ER, 4 EL, 0 Early lick
    TimeIndex[0] = 1
    Time = 1
    MaxBehave = len(behavDict['CorrectR'])
    where = Index >= MaxBehave
    if where.sum() > 0:
        logger.info('Clipped timepoints: %d, at: %s out of %d' %
                    (where.sum(), np.where(Index >= MaxBehave)[0], nRecords - 1))
    Index[Index >= MaxBehave] = MaxBehave - 1
    for i in range(1, nRecords):
        if behavDict['CorrectR'][Index[i]] and not (behavDict['EarlyLick'][Index[i]]):
            Type = 1
        elif behavDict['CorrectL'][Index[i]] and not (behavDict['EarlyLick'][Index[i]]):
            Type = 2
        elif behavDict['IncorrectR'][Index[i]] and not (behavDict['EarlyLick'][Index[i]]):
            Type = 3
        elif behavDict['IncorrectL'][Index[i]] and not (behavDict['EarlyLick'][Index[i]]):
            Type = 4
        else:
            Type = 0
        TypeIndex[i] = Type
        if Index[i] != Index[i - 1]:
            Time = 1
        else:
            Time += 1
        TimeIndex[i] = Time
    TypeIndex[0] = TypeIndex[1]
    newIndex = np.array([TypeIndex, TimeIndex]).T
    behavDict['Index'] = Index
    behavDict['newIndexAll'] = newIndex
    valid = (np.where((Index >= behavDict['start']) & (Index <= behavDict['stop'])))[0]
    behavDict['newIndex'] = newIndex[valid, :]


def aggregateByType(behavDict, data):
    """ aggregates data from an images object according to a two level index

    :param behavDict: Behavior dictionary
    :param data: Images object
    :return: a 5D numpy array (index1, index2, x, y, z)
    """

    twoLevelIndex = behavDict['newIndex']
    data.cache()
    data.count()
    si = data.shape[1:]
    maxValues = np.max(twoLevelIndex, axis=0)
    aggregated = np.zeros((maxValues[0] + 1, maxValues[1], si[0], si[1], si[2]), dtype='float32')
    aggregated[:] = np.NAN
    logger.info('aggregateByType:: ' + str(si[2]) + ' planes: ')
    # for each plane
    for k in range(0, si[2]):
        plane = data[:, :, :, k:k + 1].toarray()
        # index level 1
        for i in range(0, maxValues[0] + 1):
            # index level 2
            for j in range(0, maxValues[1]):
                idx = np.intersect1d((twoLevelIndex[:, 0] == i).nonzero()[0],
                                     (twoLevelIndex[:, 1] == (j + 1)).nonzero()[0])
                aggregated[i, j, :, :, k] = np.nanmean(plane[idx, :, :], axis=0).astype('float32')
        logger.info(k)
        sys.stdout.flush()
    return aggregated


def getEvents(volRate, behavDict):
    Events = behavDict['rawData']['saved_history']['RewardsSection_LastTrialEvents']
    Events = Events[1:]  # disregard first trial!!!
    if behavDict['start'] is not None:
        Events = Events[:behavDict['stop']]
    Events = Events[:behavDict['stop']]

    trialEvents = [None] * len(Events)

    # Get trial start frames
    stVol = np.hstack((0, np.where(np.diff(behavDict['Index']) > 0)[0]))

    def getEventStart(ID):
        Idx = np.intersect1d(np.where(evts[:, 0] == ID)[0], np.where(evts[:, 1] == 0)[0])
        if np.any(Idx):
            return evts[Idx[0], 2] - tEvents[0]
        else:
            return np.NAN

    # get the transition to state 40 from RT linux as the trigger time
    i = 0
    for trial in Events:
        tEvents = [np.NAN] * 12
        evts = np.array([np.array(xi) for xi in trial])
        triggerIdx = np.intersect1d(np.where(evts[:, 0] == 40)[0], np.where(evts[:, 1] == 0)[0])[0]
        tEvents[0] = evts[triggerIdx, 2]  # RT linux trial start time
        tEvents[1] = stVol[i]  # Imaging volume start
        if np.any(evts[:, 0] == 42):
            tEvents[2] = getEventStart(41)  # Pole Down Right
        else:
            tEvents[3] = getEventStart(41)  # Pole Down Left

        tEvents[4] = getEventStart(54)  # Start Delay
        tEvents[5] = getEventStart(55)  # Response Cue
        tEvents[6] = getEventStart(43)  # Correct Right
        tEvents[7] = getEventStart(44)  # Incorrect Right
        tEvents[8] = getEventStart(51)  # Correct Left
        tEvents[9] = getEventStart(52)  # Incorrect Left

        lickRIdx = np.where(evts[:, 1] == 3)[0]
        tEvents[10] = evts[lickRIdx, 2] - tEvents[0]  # Lick Right
        lickLIdx = np.where(evts[:, 1] == 1)[0]
        tEvents[11] = evts[lickLIdx, 2] - tEvents[0]  # Lick Left
        trialEvents[i] = tEvents
        i += 1

    behavDict['trialEventTimes'] = copy.deepcopy(trialEvents)

    for trial in trialEvents:
        for j in range(2, len(trial)):
            trial[j] = np.round(trial[j] * volRate)

    behavDict['trialEventVols'] = copy.deepcopy(trialEvents)

    for trial in trialEvents:
        for j in range(2, len(trial)):
            trial[j] = np.round(trial[j] + trial[1])

    allEventsVols = [np.NAN] * (len(trialEvents[0]) - 2)
    for trial in trialEvents:
        for j in range(len(allEventsVols)):
            allEventsVols[j] = np.hstack((allEventsVols[j], trial[j + 2]))

    allEventsVols = [(evtType[evtType < len(behavDict['Index'])]) for evtType in allEventsVols]
    allEventsVols = [(evtType[np.logical_not(np.isnan(evtType))]).astype('int') for evtType in allEventsVols]

    behavDict['allEventsVols'] = allEventsVols


def runBehavior(sc, session, regData, TForm1):
    """

    :param sc: Spark Context
    :param session: SpineSession object
    :param regData:
    :param TForm1:
    :return:
    """
    behavDict = loadBehavior(sc, session)
    aggregated = aggregateByType(behavDict, regData)
    correctR = getExampleVol(sc, aggregated[1, :140, :, :, :], TForm1, project=True)
    correctL = getExampleVol(sc, aggregated[2, :140, :, :, :], TForm1, project=True)
    session.writeTiff(correctR.transpose((1, 2, 0)), 'correctR')
    session.writeTiff(correctL.transpose((1, 2, 0)), 'correctL')
    diffs = np.diff(np.dstack(behavDict['stamps'])[0, 1, :])
    session.volRateNew = np.nanmedian(1. / diffs[diffs > 0.001])
    getEvents(session.volRateNew, behavDict)
    return behavDict
