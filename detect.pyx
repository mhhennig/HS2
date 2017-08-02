# distutils: language = c++
# distutils: sources = SpkDonline.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
import h5py
from ctypes import CDLL
import ctypes
from datetime import datetime

cdef extern from "SpkDonline.h" namespace "SpkDonline":
    cdef cppclass Detection:
        Detection() except +
        void InitDetection(long nFrames, double nSec, int sf, int NCh, long ti, long int * Indices, int agl)
        void SetInitialParams(int thres, int maa, int ahpthr, int maxsl, int minsl)
        void openSpikeFile(const char * name)
        void openFiles(const char * spikes, const char * shapes)
        void MedianVoltage(unsigned short * vm)
        void MeanVoltage(unsigned short * vm, int tInc, int tCut)
        void Iterate(unsigned short * vm, long t0, int tInc, int tCut, int tCut2)
        void FinishDetection()


def detectData(data, spikefilename, shapefilename, sfd, thres, maa = None, maxsl = None, minsl = None, ahpthr = None):
    """ Read data from a (custom, any other format would work) hdf5 file and pipe it to the spike detector. """
    # d = np.loadtxt(rawfilename)
    # f = np.load(rawfilename)
    # d = f['data']
    # d = -d + np.mean(d) + 20000
    # d = -d
    rf = h5py.File(data, 'r')
    d = rf['Raw']
    print d.shape
    if len(d.shape) != 1:
        nRecCh = d.shape[0]
    else:
        nRecCh = 1
    nFrames = d.shape[1]
    sf = int(sfd)
    nSec = nFrames / sfd  # the duration in seconds of the recording
    nSec = nFrames / sfd

    print("# Sampling rate: " + str(sf))
    print("# Number of recorded channels: " + str(nRecCh))
    print("# Analysing frames: " + str(nFrames) + ", Seconds:" +
          str(nSec))

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
    tCut = 10 + maxsl #int(0.001*int(sf)) + int(0.001*int(sf)) + 6 # what is logic behind this?
    tCut2 = 21 - maxsl
    tInc = min(nFrames-tCut-tCut2, 100000) # cap at specified number of frames
    print('tInc:'+str(tInc))
    # Messy! To be consistent, X and Y have to be swappped
    cdef np.ndarray[long, mode = "c"] Indices = np.zeros(nRecCh, dtype=ctypes.c_long)
    for i in range(nRecCh):
        Indices[i] = i
    cdef np.ndarray[unsigned short, mode="c"] vm = np.zeros((nRecCh * (tInc + tCut + tCut2)), dtype=ctypes.c_ushort)

    # initialise detection algorithm
    det.InitDetection(nFrames, nSec, sf, nRecCh, tInc, &Indices[0], int(np.mean(d)))

    det.SetInitialParams(thres, maa, ahpthr, maxsl, minsl)

    # open output file
    spikefilename = str.encode(spikefilename)
    shapefilename = str.encode(shapefilename)
    det.openFiles(spikefilename, shapefilename)

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
            vm = np.hstack((np.zeros(nRecCh * tCut), d[:,t0:t1+tCut2].flatten('F'))).astype(ctypes.c_ushort)
            print(d[:,t0:t1+tCut2].shape)
        else:
            vm = d[:,t0-tCut:t1+tCut2].flatten('F').astype(ctypes.c_ushort)
            print(d[:,t0-tCut:t1+tCut2].shape, len(vm))
        # detect spikes
        # det.MeanVoltage( &vm[0], tInc+tCut2)
        det.MeanVoltage( &vm[0], tInc, tCut)  # a bit faster (maybe)
        # print "Gets past MeanVoltage"
        det.Iterate(&vm[0], t0, tInc, tCut, tCut2)

        t0 += tInc #- tCut
        if t0 < nFrames - tCut2:
            tInc = min(tInc, nFrames - tCut2 - t0)

    det.FinishDetection()
    endTime=datetime.now()
    print('Time taken for detection: ' + str(endTime - startTime))
    print('Time per frame: ' + str(1000 * (endTime - startTime) / (nFrames)))
    print('Time per sample: ' + str(1000 *
                                    (endTime - startTime) / (nRecCh * nFrames)))
