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
        void MeanVoltage(unsigned short * vm, int tInc)
        void Iterate(unsigned short * vm, long t0, int tInc, int tCut)
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
    nFrames = 6000 #d.shape[1]
    sf = int(sfd)
    nSec = nFrames / sfd  # the duration in seconds of the recording
    tCut = 11 #int(0.001*int(sf)) + int(0.001*int(sf)) + 6 # what is logic behind this?
    nSec = nFrames / sfd
    tInc = min(nFrames-tCut, 50000) # cap at specified number of frames

    print("# Sampling rate: " + str(sf))
    print("# Number of recorded channels: " + str(nRecCh))
    print("# Analysing frames: " + str(nFrames) + ", Seconds:" +
          str(nSec))
    print("Mean: ", np.mean(d))
    print("tCut: ", tCut)

    cdef Detection * det = new Detection()

    # Messy! To be consistent, X and Y have to be swappped
    cdef np.ndarray[long, mode = "c"] Indices = np.zeros(nRecCh, dtype=ctypes.c_long)
    for i in range(nRecCh):
        Indices[i] = i
    cdef np.ndarray[unsigned short, mode="c"] vm = np.zeros((nRecCh * (tInc + tCut)), dtype=ctypes.c_ushort)

    # initialise detection algorithm
    det.InitDetection(nFrames, nSec, sf, nRecCh, tInc, &Indices[0], int(np.mean(d)))
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

    det.SetInitialParams(thres, maa, ahpthr, maxsl, minsl)

    # open output file
    spikefilename = str.encode(spikefilename)
    shapefilename = str.encode(shapefilename)
    det.openFiles(spikefilename, shapefilename)

    # cdef np.ndarray[unsigned short, ndim = 1, mode = "c"] vm =
    # np.zeros(len(d), dtype=ctypes.c_ushort)
    startTime = datetime.now()
    #vm = d.flatten('F').astype(dtype=ctypes.c_ushort)
    t0 = tCut
    while t0 + tInc <= nFrames:
        t1 = t0 + tInc
        # t1 = t0 + tInc + tCut
        print('Analysing ' + str(t1 - t0) + ' frames; ' + str(t0-tCut) + ' ' + str(t1))
        # slice data
        # this won't work when data too big to fit in memory - need different file type to .npy
        vm = d[:,t0-tCut:t1].flatten('F').astype(ctypes.c_ushort)
        print d[:,t0-tCut:t1].shape
        # detect spikes
        # det.MedianVoltage(&vm[0])
        det.MeanVoltage( &vm[0], tInc+tCut)  # a bit faster (maybe)
        det.Iterate(&vm[0], t0, tInc, tCut)

        t0 += tInc #- tCut
        if t0 < nFrames - tCut:
            tInc = min(tInc, nFrames - t0)

    det.FinishDetection()
    endTime=datetime.now()
    print('Time taken for detection: ' + str(endTime - startTime))
    print('Time per frame: ' + str(1000 * (endTime - startTime) / (nFrames)))
    print('Time per sample: ' + str(1000 *
                                    (endTime - startTime) / (nRecCh * nFrames)))
