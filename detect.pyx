# distutils: language = c++
# distutils: sources = SpkDonline.cpp SpikeHandler.cpp ProcessSpikes.cpp FilterSpikes.cpp LocalizeSpikes.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
import h5py
from ctypes import CDLL
import ctypes
from datetime import datetime
from libcpp cimport bool

cdef extern from "SpkDonline.h" namespace "SpkDonline":
    cdef cppclass Detection:
        Detection() except +
        void InitDetection(long nFrames, double nSec, int sf, int NCh, long ti, long int * Indices, int agl, int tpref, int tpostf)
        void SetInitialParams(int num_channels, int num_recording_channels, int spike_delay, int spike_peak_duration, int noise_duration, \
                              int noise_amp, int max_neighbors, int start_cutout, int end_cutout, bool to_localize, int thres, int maa, int ahpthr, int maxsl, int minsl)
        void openSpikeFile(const char * name)
        void openFiles(const char * spikes)
        void MedianVoltage(short * vm)
        void MeanVoltage(short * vm, int tInc, int tCut)
        void Iterate(short * vm, long t0, int tInc, int tCut, int tCut2)
        void FinishDetection()

def detectData(data, _num_channels, _num_recording_channels, _spike_delay, _spike_peak_duration, \
               _noise_duration, _noise_amp, _max_neighbors, _to_localize, _start_cutout, _end_cutout, \
               sfd, thres, maa = None, maxsl = None, minsl = None, ahpthr = None, tpre = 1.0, tpost = 2.2):
    """ Read data from a (custom, any other format would work) hdf5 file and pipe it to the spike detector. """
    # d = np.loadtxt(rawfilename)
    # f = np.load(rawfilename)
    # d = f['data']
    # d = -d + np.mean(d) + 20000
    # d = -d
    # rf = h5py.File(data, 'r')
    # d = rf['Raw']
    d = np.memmap(data, dtype=np.int16, mode='r')
    nRecCh = _num_channels
    nFrames = len(d)/nRecCh
    sf = int(sfd)
    nSec = nFrames / sfd  # the duration in seconds of the recording
    nSec = nFrames / sfd
    tpref = int(tpre*sf/1000)
    tpostf = int(tpost*sf/1000)
    num_channels = int(_num_channels)
    num_recording_channels = int(_num_recording_channels)
    spike_delay = int(_spike_delay)
    spike_peak_duration = int(_spike_peak_duration)
    noise_duration = int(_noise_duration)
    noise_amp = int(_noise_amp)
    max_neighbors = int(_max_neighbors)
    start_cutout = int(_start_cutout)
    end_cutout = int(_end_cutout)
    to_localize = _to_localize




    print("# Sampling rate: " + str(sf))
    print("# Number of recorded channels: " + str(nRecCh))
    print("# Analysing frames: " + str(nFrames) + ", Seconds:" +
          str(nSec))
    print("# Frames before spike in cutout: " + str(tpref))
    print("# Frames after spike in cutout: " + str(tpostf))

    cdef Detection * det = new Detection()

    # det.InitDetection(nFrames, nSec, sf, nRecCh, tInc, &Indices[0], 0)
    # set params
    # maa = None, maxsl = None, minsl = None, ahpthr = None
    if not maa:
        maa = 5
    if not maxsl:
        maxsl = int(sf*1/1000 + 0.5)
    if not minsl:
        minsl = int(sf*0.3/1000 + 0.5)
    if not ahpthr:
        ahpthr = 0

    #set tCut, tCut2 and tInc
    tCut = tpref + maxsl #int(0.001*int(sf)) + int(0.001*int(sf)) + 6 # what is logic behind this?
    tCut2 = tpostf + 1 - maxsl
    tInc = min(nFrames-tCut-tCut2, 100000) # cap at specified number of frames
    print('tInc:'+str(tInc))
    # ! To be consistent, X and Y have to be swappped
    cdef np.ndarray[long, mode = "c"] Indices = np.zeros(nRecCh, dtype=ctypes.c_long)
    for i in range(nRecCh):
        Indices[i] = i
    cdef np.ndarray[short, mode="c"] vm = np.zeros((nRecCh * (tInc + tCut + tCut2)), dtype=ctypes.c_short)
    cdef np.ndarray[short, mode = "c"] ChIndN = np.zeros((nRecCh * 10), dtype=ctypes.c_short)

    # initialise detection algorithm
    det.InitDetection(nFrames, nSec, sf, nRecCh, tInc, &Indices[0], 0, int(tpref), int(tpostf))

    det.SetInitialParams(num_channels, num_recording_channels, spike_delay, spike_peak_duration, noise_duration, \
                         noise_amp, max_neighbors, start_cutout, end_cutout, to_localize, thres, maa, ahpthr, maxsl, minsl)

    # open output file
    #spikefilename = str.encode(spikefilename)


    #det.openFiles(spikefilename)

    # cdef np.ndarray[unsigned short, ndim = 1, mode = "c"] vm =
    # np.zeros(len(d), dtype=ctypes.c_ushort)
    startTime = datetime.now()
    #vm = d.flatten('F').astype(dtype=ctypes.c_ushort)
    t0 = 0
    while t0 + tInc + tCut2 <= nFrames:
        t1 = t0 + tInc
        print('Analysing ' + str(t1 - t0) + ' frames; ' + str(t0-tCut) + ' ' + str(t1+tCut2))
        print('t0 = ' + str(t0) + ', t1 = ' +str(t1))
        # # slice data
        if t0 == 0:
            vm = np.hstack((np.zeros(nRecCh * tCut), d[:(t1+tCut2) * nRecCh])).astype(ctypes.c_short)
        else:
            vm = d[(t0-tCut) * nRecCh:(t1+tCut2) * nRecCh].astype(ctypes.c_short)
        # detect spikes
        det.MeanVoltage( &vm[0], tInc, tCut)
        det.Iterate(&vm[0], t0, tInc, tCut, tCut2)

        t0 += tInc
        if t0 < nFrames - tCut2:
            tInc = min(tInc, nFrames - tCut2 - t0)

    det.FinishDetection()
    endTime=datetime.now()
    print('Time taken for detection: ' + str(endTime - startTime))
    print('Time per frame: ' + str(1000 * (endTime - startTime) / (nFrames)))
    print('Time per sample: ' + str(1000 *
                                    (endTime - startTime) / (nRecCh * nFrames)))
