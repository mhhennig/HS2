import numpy as np
import os


def file_l(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def textfile_processor(value, s1, s2, s3, file, shp, ibc, chanlen, avg, flp):

    np.set_printoptions(threshold=np.nan)

    print("SpikeHelper: Starting loading routine...")

    squares_shape = shp
    print("SpikeHelper: ", squares_shape)

    event_bank = []
    event_num = 0
    event_sub_num = 0
    event_processor = []

    filepath = file

    file_length = file_l(filepath)

    ct = 0
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            processed_line = line.strip().split()
            if len(processed_line) != 2:
                event_sub_num += 1
            else:
                event_bank.append(int(event_sub_num))
                event_sub_num = 0
                event_num += 1
            ct += 1
            value.value = min(99.9, (100 * (float(ct) / float(file_length))))
            line = fp.readline()

    event_processor.clear()

    mmap_sizes = np.memmap(flp[0], mode='w+', shape=(3, ), dtype='int64')
    mmap_indices = np.memmap(flp[1], mode='w+', shape=(2, len(event_bank)), dtype='int32')
    mmap_data = np.memmap(flp[2], mode='w+', shape=(3, (file_length - event_num + 1)), dtype='int32')
    mmap_xy = np.memmap(flp[3], mode='w+', shape=(3, event_num), dtype='float32')
    mmap_totals = np.memmap(flp[4], mode='w+', shape=squares_shape, dtype='float32')
    mmap_spikecount = np.memmap(flp[5], mode='w+', shape=chanlen, dtype='int64')
    mmap_averages = np.memmap(flp[6], mode='w+', shape=chanlen, dtype='float64')

    s1.value = len(event_bank)
    s2.value = (file_length - event_num + 1)
    s3.value = event_num

    mmap_sizes[0] = len(event_bank)
    mmap_sizes[1] = file_length - event_num + 1
    mmap_sizes[2] = event_num

    ct = 0

    event_num = 0

    event_sub_num = 0

    with open(filepath) as fp:
        line = fp.readline()
        while line:
            processed_line = line.strip().split()
            if len(processed_line) > 2:
                mmap_data[:, event_sub_num] = np.asarray([int(processed_line[0]), int(processed_line[1]), int(processed_line[2])])
                mmap_spikecount[int(processed_line[0])] += 1
                if len(processed_line) >= 5:
                    mmap_xy[:, event_num] = np.asarray([float(processed_line[3]), float(processed_line[4]), float(processed_line[0])])
                    try:
                        if not ibc:
                            mmap_totals[int(float(processed_line[3])//10)][int(float(processed_line[4])//10)] += 1
                        else:
                            mmap_totals[int(float(processed_line[3])*2)][int(float(processed_line[4])*2)] += 1
                    except Exception:
                        pass
                event_sub_num += 1
            else:
                event_num += 1
            ct += 1
            value.value = min(99.9, (100 * (float(ct) / float(file_length))))
            line = fp.readline()

    ct = 0
    c_index = 0

    for ev in event_bank:
        ct += 1
        value.value = min(99.9, (100 * (float(ct) / float(len(event_bank)))))
        prev = c_index
        c_index += ev
        mmap_indices[:, ct - 1] = np.asarray([prev, c_index])

    event_bank.clear()

    print("SpikeHelper: Ended routine.", mmap_indices[0, 0], mmap_indices[1, 0], mmap_data[0, 0], mmap_data[1, 0], mmap_data[2, 0], mmap_xy[0, 0], mmap_xy[1, 0], len(mmap_indices[0]), len(mmap_data[0]), len(mmap_xy[0]), file_length)

    value.value = 101