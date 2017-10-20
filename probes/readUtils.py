# 3Brain 3rd gen .brw (HDF5)
import h5py
import ctypes


def read_flat(d, t0, t1, nch):
    return d[t0*nch:t1*nch].astype(ctypes.c_short)


def openHDF5file(path):
    return h5py.File(path, 'r')


def getHDF5params(rf):
    # Read recording variables
    recVars = rf.require_group('3BRecInfo/3BRecVars/')
    # bitDepth = recVars['BitDepth'].value[0]
    # maxV = recVars['MaxVolt'].value[0]
    # minV = recVars['MinVolt'].value[0]
    nFrames = recVars['NRecFrames'].value[0]
    samplingRate = recVars['SamplingRate'].value[0]
    # signalInv = recVars['SignalInversion'].value[0]

    # Read chip variables
    chipVars = rf.require_group('3BRecInfo/3BMeaChip/')
    # nRows = chipVars['NRows'].value[0]
    nCols = chipVars['NCols'].value[0]
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

    print('3Brain data format: '+str(file_format))
    # Compute indices
    rawIndices = rf['3BRecInfo/3BMeaStreams/Raw/Chs'].value

    # Name channels ([0..4095] for fullarray files)
    chIndices = [(x-1) + (y-1)*nCols for (y, x) in rawIndices]
    # chIndices = [(x-1) + (y-1)*nCols for (x,y) in rawIndices]
    # Swap X and Y (old format)

    return (nFrames, samplingRate, nRecCh, chIndices, file_format)


def readHDF5(rf, t0, t1):
    ''' In order to use the algorithms designed for the old format,
    the input data must be inverted.'''
    return 4095 - rf['3BData/Raw'][t0:t1].flatten().astype(ctypes.c_short)


def readHDF5t_100(rf, t0, t1, nch):
    ''' Transposed version for the interpolation method. '''
    if t0 <= t1:
        d = 2048 - rf['3BData/Raw'][t0:t1].flatten('C').astype(ctypes.c_short)
        d[d > 200] = 0
        return d
    else:  # Reversed read
        raise Exception('Not sure about this one.')

        return 2048 - rf['3BData/Raw'][t1:t0].flatten(
                    'F').astype(ctypes.c_short)


def readHDF5t_101(rf, t0, t1, nch):
    ''' Transposed version for the interpolation method. '''
    if t0 <= t1:
        return rf['3BData/Raw'][nch*t0:nch*t1].reshape(
            (-1, nch), order='C').flatten('F').astype(ctypes.c_short)
    else:  # Reversed read
        return rf['3BData/Raw'][nch*t1:nch*t0].reshape(
            (-1, nch), order='C').flatten('F').astype(ctypes.c_short)
