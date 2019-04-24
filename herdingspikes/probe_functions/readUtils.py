# 3Brain 3rd gen .brw (HDF5)
import h5py
import ctypes

import numpy as np


def read_flat(d, t0, t1, nch):
    return d[t0*nch:t1*nch].astype(ctypes.c_short)


def openHDF5file(path, **kwargs):
    return h5py.File(path, 'r', **kwargs)


def getHDF5params(rf):
    # Read recording variables
    recVars = rf.require_group('3BRecInfo/3BRecVars/')
    # bitDepth = recVars['BitDepth'].value[0]
    # maxV = recVars['MaxVolt'].value[0]
    # minV = recVars['MinVolt'].value[0]
    nFrames = recVars['NRecFrames'][()][0]
    samplingRate = recVars['SamplingRate'][()][0]
    signalInv = recVars['SignalInversion'][()][0]

    # Read chip variables
    chipVars = rf.require_group('3BRecInfo/3BMeaChip/')
    # nRows = chipVars['NRows'].value[0]
    nCols = chipVars['NCols'][()][0]
    # nChipCh = nRows * nCols # Total number of channels

    # Get the actual number of channels used in the recording
    file_format = rf['3BData'].attrs.get('Version')
    if file_format == 100:
        nRecCh = len(rf['3BData/Raw'][0])
        # raise Warning('This may go wrong!')
    elif file_format == 101:
        nRecCh = int(1.*rf['3BData/Raw'].shape[0]/nFrames)
    else:
        raise Exception('Unknown data file format.')

    print('# 3Brain data format:', file_format, 'signal inversion', signalInv)
    print('#       signal range: ', recVars['MinVolt'][()][0], '- ',
          recVars['MaxVolt'][()][0])
    # Compute indices
    rawIndices = rf['3BRecInfo/3BMeaStreams/Raw/Chs'][()]

    # Name channels ([0..4095] for fullarray files)
    chIndices = [(x-1) + (y-1)*nCols for (y, x) in rawIndices]
    # chIndices = [(x-1) + (y-1)*nCols for (x,y) in rawIndices]
    # Swap X and Y (old format)

    return (nFrames, samplingRate, nRecCh, chIndices, file_format, signalInv)


def readHDF5(rf, t0, t1):
    ''' In order to use the algorithms designed for the old format,
    the input data must be inverted.'''
    return 4095 - rf['3BData/Raw'][t0:t1].flatten().astype(ctypes.c_short)


def readHDF5t_100(rf, t0, t1, nch):
    ''' Transposed version for the interpolation method. '''
    if t0 <= t1:
        d = 2048 - rf['3BData/Raw'][t0:t1].flatten('C').astype(ctypes.c_short)
        d[np.where(np.abs(d) > 1500)[0]] = 0
        return d
    else:  # Reversed read
        raise Exception('Reading backwards? Not sure about this.')
        return 2048 - rf['3BData/Raw'][t1:t0].flatten(
                    'F').astype(ctypes.c_short)


def readHDF5t_100_i(rf, t0, t1, nch):
    ''' Transposed version for the interpolation method. '''
    if t0 <= t1:
        d = rf['3BData/Raw'][t0:t1].flatten('C').astype(ctypes.c_short) - 2048
        d[np.where(np.abs(d) > 1500)[0]] = 0
        return d
    else:  # Reversed read
        raise Exception('Reading backwards? Not sure about this.')
        return rf['3BData/Raw'][t1:t0].flatten(
                    'F').astype(ctypes.c_short) - 2048


def readHDF5t_101(rf, t0, t1, nch):
    ''' Transposed version for the interpolation method. '''
    if t0 <= t1:
        d = rf['3BData/Raw'][nch*t0:nch*t1].reshape(
            (-1, nch), order='C').flatten('C').astype(ctypes.c_short)-2048
        d[np.abs(d) > 1500] = 0
        return d
    else:  # Reversed read
        raise Exception('Reading backwards? Not sure about this.')
        d = rf['3BData/Raw'][nch*t1:nch*t0].reshape(
            (-1, nch), order='C').flatten('C').astype(ctypes.c_short)-2048
        d[np.where(np.abs(d) > 1500)[0]] = 0
        return d


def readHDF5t_101_i(rf, t0, t1, nch):
    ''' Transposed version for the interpolation method. '''
    if t0 <= t1:
        d = 2048-rf['3BData/Raw'][nch*t0:nch*t1].reshape(
            (-1, nch), order='C').flatten('C').astype(ctypes.c_short)
        d[np.where(np.abs(d) > 1500)[0]] = 0
        return d
    else:  # Reversed read
        raise Exception('Reading backwards? Not sure about this.')
        d = 2048-rf['3BData/Raw'][nch*t1:nch*t0].reshape(
            (-1, nch), order='C').flatten('C').astype(ctypes.c_short)
        d[np.where(np.abs(d) > 1500)[0]] = 0
        return d


def getNeuroSeekerParams(rf, pipette=False):
    if pipette:
        return rf['kampff_probe_data'].attrs['no_frames'], \
            rf.attrs['frequency'], 1, [1]

    return (rf['kampff_probe_data'].attrs['no_frames'],
            rf.attrs['frequency'],
            rf['kampff_probe_data'].attrs['no_channels'],
            np.arange(rf['kampff_probe_data'].attrs['no_channels']),
            rf.attrs['probe_closest_electrode'])


def readNeuroSeekerProbe(rf, t0, t1):
    return (rf['kampff_probe_data'][t0:t1].flatten() - 20000).astype(
        ctypes.c_short)


def readNeuroSeekerPipette(rf, t0, t1):
    return 50000 - rf['kampff_pipette_data'][t0:t1].flatten()


def readSiNAPS_S1Probe(raw_data, t0, t1):
    raw_traces = raw_data[t0:t1]
    return raw_traces.flatten('C').astype(ctypes.c_short)
