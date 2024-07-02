# distutils: language = c++

import cython
import numpy as np
cimport numpy as np
cimport cython
import h5py
from ctypes import CDLL
import ctypes
from datetime import datetime
from libcpp cimport bool
from libcpp.string cimport string
import sys
import pprint
import warnings

cdef extern from "SpkDonline.h" namespace "SpkDonline":
    cdef cppclass Detection:
        Detection() except +
        void InitDetection(long nFrames, int sf, int NCh, long ti, long int * Indices, int agl)
        void SetInitialParams(float * pos_mtx, int * neigh_mtx, int num_channels,
                              int spike_peak_duration, string file_name, int noise_duration,
                              float noise_amp_percent, float inner_radius, int* _masked_channels, \
                              int max_neighbors, int num_com_centers, bool to_localize, int thres, int cutout_start, int cutout_end, \
                              int maa, int ahpthr, int maxsl, int minsl, bool decay_filtering, bool verbose)
        void MedianVoltage(short * vm, int tInc, int tCut)
        void MeanVoltage(short * vm, int tInc, int tCut)
        void Iterate(short *vm, long t0, int tInc, int tCut, int tCut2, int maxFramesProcessed)
        void FinishDetection()

def detectData(probe, file_name, to_localize, sf, thres,
               cutout_start, cutout_end,
               maa=5, maxsl=None, minsl=None,
               ahpthr=0, num_com_centers=1,
               decay_filtering=False, verbose=True,
               nFrames=None, tInc=50000):
    """ Read data from a file and pipe it to the spike detector. """

    # READ PROBE PARAMETERS
    sf = int(sf) # ensure sampling rate is integer, assumed to be in Hertz
    num_channels = int(probe.num_channels)
    spike_peak_duration = int(probe.spike_peak_duration)
    noise_duration = int(probe.noise_duration)
    noise_amp_percent = float(probe.noise_amp_percent)
    max_neighbors = int(probe.max_neighbors)
    inner_radius = float(probe.inner_radius)
    # positions_file_path = probe.positions_file_path.encode()
    # neighbors_file_path = probe.neighbors_file_path.encode()

    if nFrames is None:
        nFrames = probe.nFrames

    # READ DETECTION PARAMETERS AND SET DEFAULTS
    nRecCh = num_channels
    if not maxsl:
        maxsl = int(sf*1/1000 + 0.5)
    if not minsl:
        minsl = int(sf*0.3/1000 + 0.5)

    masked_channel_list = probe.masked_channels
    cdef np.ndarray[int, mode="c"] masked_channels = np.ones(num_channels, dtype=ctypes.c_int)
    if masked_channel_list == []:
        print("# Not Masking any Channels")
        masked_channel_list = None
    if masked_channel_list is not None:
        print("# Masking Channels: " +str(masked_channel_list))
        for channel in masked_channel_list:
            masked_channels[channel] = 0

    print("# Sampling rate: " + str(sf))
    if to_localize:
        print("# Localization On")
    else:
        print("# Localization Off")
    if verbose:
        print("# Writing out extended detection info")
    print("# Number of recorded channels: " + str(num_channels))
    if num_channels<20:
        print("# Few recording channels: not subtracing mean from activity")
    print("# Analysing frames: " + str(nFrames) + "; Seconds: " + str(nFrames/sf))
    print("# Frames before spike in cutout: " + str(cutout_start))
    print("# Frames after spike in cutout: " + str(cutout_end))

    cdef Detection * det = new Detection()

    # set tCut, tCut2 and tInc
    tCut = cutout_start + maxsl
    tCut2 = cutout_end + maxsl
    print("# tcuts: " + str(tCut) + " "+ str(tCut2) )

    tInc = min(nFrames-tCut-tCut2, tInc) # cap at specified number of frames
    maxFramesProcessed = tInc;
    print('# tInc: '+str(tInc))
    # ! To be consistent, X and Y have to be swappped
    cdef np.ndarray[long, mode = "c"] Indices = np.zeros(nRecCh, dtype=ctypes.c_long)
    for i in range(nRecCh):
        Indices[i] = i
    cdef np.ndarray[short, mode = "c"] vm = np.zeros((nRecCh * (tInc + tCut + tCut2)), dtype=ctypes.c_short)
    #cdef np.ndarray[short, mode = "c"] ChIndN = np.zeros((nRecCh * 10), dtype=ctypes.c_short)

    # initialise detection algorithm
    det.InitDetection(nFrames, sf, nRecCh, tInc, &Indices[0], 0)

    cdef np.ndarray[float, ndim=2, mode = "c"] position_matrix = np.zeros((nRecCh,2), dtype=ctypes.c_float)
    for i,p in enumerate(probe.positions):
      position_matrix[i,0] = p[0]
      position_matrix[i,1] = p[1]
    cdef np.ndarray[int, ndim=2, mode = "c"] neighbor_matrix = np.zeros((
        nRecCh,np.max([len(p) for p in probe.neighbors])), dtype=ctypes.c_int)-1
    for i,p in enumerate(probe.neighbors):
      neighbor_matrix[i,:len(p)] = p

    det.SetInitialParams(&position_matrix[0,0], &neighbor_matrix[0,0], num_channels,
                         spike_peak_duration, file_name, noise_duration,
                         noise_amp_percent, inner_radius, &masked_channels[0],
                         max_neighbors, num_com_centers, to_localize,
                         thres, cutout_start, cutout_end, maa, ahpthr, maxsl,
                         minsl, decay_filtering, verbose)

    startTime = datetime.now()
    t0 = 0
    while t0 + tInc + tCut2 <= nFrames:
        t1 = t0 + tInc
        if verbose:
            print('# Analysing frames from ' + str(t0-tCut) + ' to '+str(t1+tCut2)+\
              '  ({:.1f}%)'.format(100*t0/nFrames))
        #sys.stdout.flush()
        # slice data
        if t0 == 0:
            vm = np.hstack((np.zeros(nRecCh * tCut, dtype=ctypes.c_short), probe.Read(0, t1+tCut2))).astype(ctypes.c_short)
        else:
            vm = probe.Read(t0-tCut, t1+tCut2)
        # detect spikes
        # if num_channels>=20:
        #     det.MeanVoltage( &vm[0], tInc, tCut)
        det.Iterate(&vm[0], t0, tInc, tCut, tCut2, maxFramesProcessed)
        t0 += tInc
        if t0 < nFrames - tCut2:
            tInc = min(tInc, nFrames - tCut2 - t0)

    now = datetime.now()
    #Save state of detection
    detection_state_dict = {
            'Probe Name': probe.__class__.__name__,
            'Probe Object': probe,
            'Date and Time Detection': str(now),
            'Threshold': thres,
            'Localization': to_localize,
            'Masked Channels': masked_channel_list,
            'Associated Results File': file_name,
            # 'Positions File Path': positions_file_path,
            # 'Neighbors File Path': neighbors_file_path,
            'Cutout Length': cutout_start + cutout_end,
            'Advice': 'For more information about detection, load and look at the parameters of the probe object',
        }

    target = open(file_name.decode() + 'DetectionDict' + now.strftime("%Y-%m-%d_%H%M%S_%f") + '.txt', 'a')
    target.write(pprint.pformat(detection_state_dict))
    target.close()

    det.FinishDetection()
    endTime=datetime.now()
    print('# Detection completed, time taken: ' + str(endTime - startTime))
    print('# Time per frame: ' + str(1000 * (endTime - startTime) / (nFrames)))
    print('# Time per sample: ' + str(1000 * (endTime - startTime) / (nRecCh * nFrames)))
