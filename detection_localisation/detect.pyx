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
from libcpp.string cimport string
import sys

cdef extern from "SpkDonline.h" namespace "SpkDonline":
    cdef cppclass Detection:
        Detection() except +
        void InitDetection(long nFrames, double nSec, int sf, int NCh, long ti, long int * Indices, int agl, int tpref, int tpostf)
        void SetInitialParams(string positions_file_path, string neighbors_file_path, int num_channels, int num_recording_channels, int spike_delay,
                              int spike_peak_duration, string file_name, int noise_duration, float noise_amp_percent, int* _masked_channels, \
                              int max_neighbors, int inner_radius, bool to_localize, int thres, int cutout_start, int cutout_end, \
                              int maa, int ahpthr, int maxsl, int minsl)
        void MedianVoltage(short * vm)
        void MeanVoltage(short * vm, int tInc, int tCut)
        void Iterate(short * vm, long t0, int tInc, int tCut, int tCut2, int maxFramesProcessed)
        void FinishDetection()


def read_flat(d, t0, t1, nch):
  return d[t0*nch:t1*nch].astype(ctypes.c_short)


def detectData(probe, _file_name, _to_localize, sfd, thres,
               _cutout_start=10, _cutout_end=20, maa=5, maxsl=None, minsl=None,
               ahpthr=0, tpre=1.0, tpost=2.2):
    """ Read data from a file and pipe it to the spike detector. """

    # if data_format is 'flat':
    #   d = np.memmap(filename, dtype=np.int16, mode='r')
    #   nRecCh = _num_channels
    #   assert len(d)/nRecCh==len(d)//nRecCh, 'data not multiple of channel number'
    #   nFrames = len(d)//nRecCh
    #   sf = int(sfd)
    #   read_function = read_flat
    # elif data_format is 'biocam':
    #   from readUtils import openHDF5file, getHDF5params, readHDF5t_100, readHDF5t_101
    #   d = openHDF5file(filename)
    #   nFrames, sfd, nRecCh, chIndices, file_format = getHDF5params(d)
    #   _num_channels = nRecCh
    #   _num_recording_channels = nRecCh
    #
    #   if file_format == 100:
    #       read_function = readHDF5t_100
    #   else:
    #       read_function = readHDF5t_101
    # else:
    #   raise NotImplementedError('Unknown file format')

    # nFrames = 14000

    nSec = probe.nFrames / sfd  # the duration in seconds of the recording
    nSec = probe.nFrames / sfd
    sf = int(sfd)
    tpref = int(tpre*sf/1000)
    tpostf = int(tpost*sf/1000)
    num_channels = int(probe.num_channels)
    num_recording_channels = int(probe.num_recording_channels)
    spike_delay = int(probe.spike_delay)
    spike_peak_duration = int(probe.spike_peak_duration)
    noise_duration = int(probe.noise_duration)
    noise_amp_percent = float(probe.noise_amp_percent)
    max_neighbors = int(probe.max_neighbors)
    inner_radius = int(probe.inner_radius)
    cutout_start = int(_cutout_start)
    cutout_end = int(_cutout_end)
    to_localize = _to_localize
    nRecCh = num_channels
    nFrames = probe.nFrames
    masked_channel_list = probe.masked_channels
    cdef np.ndarray[int, mode="c"] masked_channels = np.ones(num_channels, dtype=ctypes.c_int)
    if masked_channel_list == []:
        masked_channel_list = None
    if masked_channel_list != None:
        for channel in masked_channel_list:
            masked_channels[channel] = 0
    for channel in range(num_recording_channels, num_channels):
        masked_channels[channel] = 0

    # positions_file_path = str(_positions_file_path)
    # neighbors_file_path = str(_neighbors_file_path)
    positions_file_path = probe.positions_file_path.encode() # <- python 3 seems to need this
    neighbors_file_path = probe.neighbors_file_path.encode()

    print("# Sampling rate: " + str(sf))

    if to_localize == True:
        print("# Localization On")
    else:
        print("# Localization Off")

    if masked_channel_list != None:
        print("# Masking Channels: " +str(masked_channel_list))
    else:
        print("# Not Masking any Channels")

    print("# Number of recorded channels: " + str(num_channels))
    print("# Analysing frames: " + str(nFrames) + ", Seconds:" +
          str(nSec))
    print("# Frames before spike in cutout: " + str(tpref))
    print("# Frames after spike in cutout: " + str(tpostf))

    cdef Detection * det = new Detection()

    if not maxsl:
        maxsl = int(sf*1/1000 + 0.5)
    if not minsl:
        minsl = int(sf*0.3/1000 + 0.5)

    # set tCut, tCut2 and tInc
    # tCut = tpref + maxsl #int(0.001*int(sf)) + int(0.001*int(sf)) + 6 # what is logic behind this?
    # FIX: make sure enough data points are available:
    tCut = max((tpref + maxsl, cutout_start+maxsl))
    # tCut2 = tpostf + 1 - maxsl
    # FIX: make sure enough data points are available:
    tCut2 = max(tpostf + 1 - maxsl, cutout_end+maxsl)
    print("# tcuts: " + str(tCut) + " "+ str(tCut2) )

    tInc = min(nFrames-tCut-tCut2, 200000) # cap at specified number of frames
    # tInc = 10000
    maxFramesProcessed = tInc;
    print('# tInc: '+str(tInc))
    # ! To be consistent, X and Y have to be swappped
    cdef np.ndarray[long, mode = "c"] Indices = np.zeros(nRecCh, dtype=ctypes.c_long)
    for i in range(nRecCh):
        Indices[i] = i
    cdef np.ndarray[short, mode="c"] vm = np.zeros((nRecCh * (tInc + tCut + tCut2)), dtype=ctypes.c_short)
    cdef np.ndarray[short, mode = "c"] ChIndN = np.zeros((nRecCh * 10), dtype=ctypes.c_short)

    # initialise detection algorithm
    det.InitDetection(nFrames, nSec, sf, nRecCh, tInc, &Indices[0], 0, int(tpref), int(tpostf))

    det.SetInitialParams(positions_file_path, neighbors_file_path, num_channels, num_recording_channels, spike_delay, spike_peak_duration,
                         _file_name, noise_duration, noise_amp_percent, &masked_channels[0], max_neighbors, inner_radius, to_localize, thres,
                         cutout_start, cutout_end, maa, ahpthr, maxsl, minsl)
    startTime = datetime.now()
    t0 = 0
    while t0 + tInc + tCut2 <= nFrames:
        t1 = t0 + tInc
        print('# Analysing ' + str(t1 - t0) + ' frames; ' + str(t0-tCut) + ' ' + str(t1+tCut2))
        # print('# t0 = ' + str(t0) + ', t1 = ' +str(t1)+', from '+str(t0-tCut))
        sys.stdout.flush()
        # # slice data
        if t0 == 0:
            #vm = np.hstack((np.zeros(nRecCh * tCut), d[:(t1+tCut2) * nRecCh])).astype(ctypes.c_short)
            vm = np.hstack((np.zeros(nRecCh * tCut, dtype=ctypes.c_short), probe.Read(0, t1+tCut2)))
        else:
            #vm = d[(t0-tCut) * nRecCh:(t1+tCut2) * nRecCh].astype(ctypes.c_short)
            vm = probe.Read(t0-tCut, t1+tCut2)
        # detect spikes
        det.MeanVoltage( &vm[0], tInc, tCut)
        det.Iterate(&vm[0], t0, tInc, tCut, tCut2, maxFramesProcessed)
        t0 += tInc
        if t0 < nFrames - tCut2:
            tInc = min(tInc, nFrames - tCut2 - t0)

    det.FinishDetection()
    endTime=datetime.now()
    print('# Time taken for detection: ' + str(endTime - startTime))
    print('# Time per frame: ' + str(1000 * (endTime - startTime) / (nFrames)))
    print('# Time per sample: ' + str(1000 *
(endTime - startTime) / (nRecCh * nFrames)))
