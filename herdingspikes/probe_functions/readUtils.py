# Various helper functions to read from BioCam data files
# 3Brain 3rd gen .brw (HDF5)

import h5py
import ctypes
import json
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
    elif (file_format == 101) or (file_format == 102):
        nRecCh = int(1. * rf['3BData/Raw'].shape[0] / nFrames)
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

def getHDF5params_brw4(rf):

    exp_setting = json.loads(rf['ExperimentSettings'].asstr()[0])
    
    samplingRate = exp_setting['TimeConverter']['FrameRate']
    signalInv = exp_setting['ValueConverter']['ScaleFactor']
    file_format = 'brw4'

    for key in rf:
        if key.startswith("Well_"):
            nFrames = int((len(rf[key]['Raw']))/len(rf[key]['StoredChIdxs']))
            # Get the actual number of channels used in the recording
            nRecCh = len(rf[key]['StoredChIdxs'])
            # Compute indices
            chIndices = rf[key]['StoredChIdxs'][()].tolist()

    print('# 3Brain data format:', 'BRW v4.x', 'signal inversion', signalInv)
    print('#       signal range: ', exp_setting['ValueConverter']['MinAnalogValue'], '- ',
          exp_setting['ValueConverter']['MaxAnalogValue'])
    
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

def readHDF5_brw4(rf, t0, t1, nch):
    ''' Transposed version for the interpolation method. '''
    for key in rf:
        if key.startswith("Well_"):
            if t0 <= t1:
                d = rf[key]['Raw'][nch*t0:nch*t1].reshape(
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
