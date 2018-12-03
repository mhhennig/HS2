from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from .probe_functions.readUtils import read_flat, readSiNAPS_S1Probe
from .probe_functions.readUtils import openHDF5file, getHDF5params
from .probe_functions.readUtils import readHDF5t_100, readHDF5t_101
from .probe_functions.readUtils import readHDF5t_100_i, readHDF5t_101_i
from .probe_functions.readUtils import getNeuroSeekerParams, readNeuroSeekerProbe
from .probe_functions.neighborMatrixUtils import createNeighborMatrix
import h5py
import ctypes
import os.path

def in_probes_dir(file):
    return os.path.join(this_file, 'probes', file)

from scipy.spatial.distance import cdist

def create_probe_files(pos_file, neighbor_file, radius, ch_positions):
    n_channels = ch_positions.shape[0]
    # NB: Notice the column, row order in write
    with open(pos_file, 'w') as f:
        for pos in ch_positions:
            f.write('{},{},\n'.format(int(pos[0]), int(pos[1])))
    f.close()
    # # NB: it is also possible to use metric='cityblock' (Manhattan distance)
    distances = cdist(ch_positions, ch_positions, metric='euclidean')
    indices = np.arange(n_channels)
    with open(neighbor_file, "w") as f:
        for dist_from_ch in distances:
            neighbors = indices[dist_from_ch <= radius]
            f.write('{},\n'.format(str(list(neighbors))[1:-1]))
    f.close()

this_file_path = os.path.dirname(os.path.abspath(__file__))
def in_probe_info_dir(file):
    return os.path.join(this_file_path, 'probe_info', file)

class NeuralProbe(object):

    def __init__(self, num_channels, spike_delay, spike_peak_duration,
                 noise_duration, noise_amp_percent, inner_radius,
                 fps, positions_file_path, neighbors_file_path,
                 neighbor_radius, masked_channels,
                 recorded_channels=None):
        if(neighbor_radius is not None):
            createNeighborMatrix(neighbors_file_path, positions_file_path, neighbor_radius)
        self.num_channels = num_channels
        self.spike_delay = spike_delay
        self.spike_peak_duration = spike_peak_duration
        self.noise_duration = noise_duration
        self.noise_amp_percent = noise_amp_percent
        self.fps = fps
        self.positions_file_path = positions_file_path
        self.neighbors_file_path = neighbors_file_path
        self.masked_channels = masked_channels
        self.inner_radius = inner_radius
        if masked_channels is None:
            self.masked_channels = []

        self.loadPositions(positions_file_path)
        self.loadNeighbors(neighbors_file_path)

        # if a probe only records a subset of channels, filter out unused channels
        # this may happen in Biocam recordings
        # note this uses positions to identify recorded channels
        # requires channels are an ordered list
        if recorded_channels is not None:
            inds = np.zeros(self.num_channels, dtype=int)
            # recorded_channels = np.roll(recorded_channels,1,axis=1)
            for i, c in enumerate(recorded_channels):
                inds[i] = np.where(np.all((
                    self.positions-np.array([c[1]-1, c[0]-1])) == 0, axis=1))[0]
            self.positions = self.positions[inds]
            x0 = np.min([p[0] for p in self.positions])
            y0 = np.min([p[1] for p in self.positions])
            x1 = np.max([p[0] for p in self.positions])
            y1 = np.max([p[1] for p in self.positions])
            print('# ', x0, x1, y0, y1, x1-x0+1, y1-y0+1, num_channels)
            lm = np.zeros((64, 64), dtype=int)-1
            # oddness because x/y are transposed in brw
            lm[y0:y1+1, x0:x1+1] = np.arange(self.num_channels).T.reshape(
                y1-y0+1, x1-x0+1)
            self.neighbors = [lm.flatten()[self.neighbors[i]] for i in inds]
            self.neighbors = [n[(n >= 0)] for n in self.neighbors]
            self.positions = self.positions-np.min(self.positions, axis=0)

    # Load in neighbor and positions files
    def loadNeighbors(self, neighbors_file_path):
        neighbor_file = open(neighbors_file_path, 'r')
        neighbors = []
        for neighbor in neighbor_file.readlines():
            neighbors.append(np.array(neighbor[:-2].split(',')).astype(int))
        neighbor_file.close()
        # assert len(neighbors) == len(pos)
        self.neighbors = neighbors
        self.max_neighbors = max([len(n) for n in neighbors])

    def loadPositions(self, positions_file_path):
        position_file = open(positions_file_path, 'r')
        positions = []
        for position in position_file.readlines():
            positions.append(np.array(position[:-2].split(',')).astype(int))
        self.positions = np.asarray(positions)
        position_file.close()

    # Show visualization of probe
    def show(self, show_neighbors=[10], figwidth=3):
        xmax, ymax = self.positions.max(0)
        xmin, ymin = self.positions.min(0)
        ratio = ymax / xmax
        plt.figure(figsize=(figwidth, figwidth * ratio))
        for ch in show_neighbors:
            for neighbor in self.neighbors[ch]:
                plt.plot([self.positions[ch, 0], self.positions[neighbor, 0]],
                         [self.positions[ch, 1], self.positions[neighbor, 1]],
                         '--k', alpha=.7)
        plt.scatter(*self.positions.T)
        plt.scatter(*self.positions[self.masked_channels].T, c='r')
        for i, pos in enumerate(self.positions):
            plt.annotate(i, pos)

    def Read(self, t0, t1):
        raise NotImplementedError("The Read function is not implemented for \
            this probe")

    def getChannelsPositions(self, channels):
        channel_positions = []
        for channel in channels:
            if channel >= self.num_channels:
                raise ValueError(
                    'Channel index too large, maximum ' + self.num_channels)
            else:
                channel_positions.append(self.positions[channel])
        return channel_positions


class NeuroPixel(NeuralProbe):
    def __init__(self, data_file_path=None, fps=30000, num_channels=385,
                 spike_delay=5, spike_peak_duration=4, noise_duration=3,
                 noise_amp_percent=1, inner_radius=76, neighbor_radius=None,
                 masked_channels=None):
        positions_file_path = in_probe_info_dir('positions_neuropixel')
        neighbors_file_path = in_probe_info_dir('neighbormatrix_neuropixel')
        NeuralProbe.__init__(
            self,
            num_channels=num_channels,
            spike_delay=spike_delay,
            spike_peak_duration=spike_peak_duration,
            noise_duration=noise_duration,
            noise_amp_percent=noise_amp_percent,
            fps=fps,
            inner_radius=inner_radius,
            positions_file_path=positions_file_path,
            neighbors_file_path=neighbors_file_path,
            neighbor_radius=neighbor_radius,
            masked_channels=masked_channels
            )

        self.data_file = data_file_path
        if data_file_path is not None:
            self.d = np.memmap(data_file_path, dtype=np.int16, mode='r')
            assert len(self.d) / self.num_channels == len(self.d) // \
                self.num_channels, 'Data not multiple of channel number'
            self.nFrames = len(self.d) // self.num_channels
        else:
            print('Note: data file not specified, things may break')

    def Read(self, t0, t1):
        return read_flat(self.d, t0, t1, self.num_channels)

    def show(self, show_neighbors=[10], figwidth=3):
        masked_channels = self.masked_channels[:-1]
        positions = self.positions[:-1]
        xmax, ymax = positions.max(0)
        xmin, ymin = positions.min(0)
        ratio = ymax / xmax
        plt.figure(figsize=(figwidth, figwidth * ratio))
        for ch in show_neighbors:
            for neighbor in self.neighbors[ch]:
                plt.plot([positions[ch, 0], positions[neighbor, 0]],
                         [positions[ch, 1], positions[neighbor, 1]],
                         '--k', alpha=.7)
        plt.scatter(*positions.T)
        plt.scatter(*positions[masked_channels].T, c='r')
        for i, pos in enumerate(positions):
            plt.annotate(i, pos)


class BioCam(NeuralProbe):
    def __init__(self, data_file_path=None, fps=0, spike_delay=5,
                 spike_peak_duration=4, noise_duration=4, noise_amp_percent=1,
                 inner_radius=1.75, neighbor_radius=None, masked_channels=[0]):
        self.data_file = data_file_path
        if data_file_path is not None:
            self.d = openHDF5file(data_file_path)
            self.nFrames, sfd, nRecCh, chIndices, file_format, inversion = getHDF5params(
                self.d)
            print("# Signal inversion looks like", str(inversion), ", guessing the "
                  "right method for data access.\n# If your detection results "
                  "look strange, signal polarity is wrong.\n")
            if file_format == 100:
                if inversion == -1:
                    self.read_function = readHDF5t_100
                else:
                    self.read_function = readHDF5t_100_i
            else:
                if inversion == -1:
                    self.read_function = readHDF5t_101_i
                else:
                    self.read_function = readHDF5t_101
        else:
            print('# Note: data file not specified, setting some defaults')
            nRecCh = 4096
            sfd = fps
        if nRecCh < 4096:
            print('# Note: only', nRecCh, 'channels recorded, fixing positions/neighbors')
            print('# This may break - known to work only for rectangular sections!')
            recorded_channels = self.d['3BRecInfo']['3BMeaStreams']['Raw']['Chs']
        else:
            recorded_channels = None

        positions_file_path = in_probe_info_dir('positions_biocam')
        neighbors_file_path = in_probe_info_dir('neighbormatrix_biocam')
        NeuralProbe.__init__(
            self,
            num_channels=nRecCh,
            spike_delay=spike_delay,
            spike_peak_duration=spike_peak_duration,
            noise_duration=noise_duration,
            noise_amp_percent=noise_amp_percent,
            fps=sfd,
            inner_radius=inner_radius,#1.75,
            positions_file_path=positions_file_path,
            neighbors_file_path=neighbors_file_path,
            neighbor_radius=neighbor_radius,
            masked_channels=masked_channels,
            recorded_channels=recorded_channels
            )
    def Read(self, t0, t1):
        return self.read_function(self.d, t0, t1, self.num_channels)

class MCS120(NeuralProbe):
    def __init__(self, data_file_path=None, fps=10000, spike_delay=5,
                 spike_peak_duration=4, noise_amp_percent=1, inner_radius=1.75,
                 neighbor_radius=None, masked_channels=None):
        self.data_file = data_file_path
        print("# BETA: Initialising MCS120 probe. Note nothing is known about the geometry. ")
        if data_file_path is not None:
            d = h5py.File(data_file_path)
            self.d = d
            nRecCh = d['Data']['Recording_0']['AnalogStream']['Stream_0']['ChannelData'].shape[0]
            self.nFrames = d['Data']['Recording_0']['AnalogStream']['Stream_0']['ChannelData'].shape[1]
            sfd = fps
        else:
            print('# Note: data file not specified, setting some defaults')
            nRecCh = 120
            sfd = fps

        positions_file_path = in_probe_info_dir('positions_mcs120')
        neighbors_file_path = in_probe_info_dir('neighbormatrix_mcs120')
        NeuralProbe.__init__(
            self,
            num_channels=nRecCh,
            spike_delay=spike_delay,
            spike_peak_duration=spike_peak_duration,
            noise_duration=int(sfd*0.0009),
            noise_amp_percent=noise_amp_percent,
            fps=sfd,
            inner_radius=inner_radius,
            positions_file_path=positions_file_path,
            neighbors_file_path=neighbors_file_path,
            neighbor_radius=neighbor_radius,
            masked_channels=masked_channels
            )

    def Read(self, t0, t1):
        return self.d['Data']['Recording_0']['AnalogStream']['Stream_0']['ChannelData'][:,t0:t1].T.ravel().astype(ctypes.c_short)


class Mea1k(NeuralProbe):
    def __init__(self, data_file_path=None, fps=20000, number_of_frames=4450600,
                 num_channels=69, spike_delay=5, spike_peak_duration=4,
                 noise_duration=2, noise_amp_percent=1, inner_radius=80,
                 neighbor_radius=None, masked_channels=None):
        positions_file_path = in_probe_info_dir('positions_mea1k')
        neighbors_file_path = in_probe_info_dir('neighbormatrix_mea1k')
        NeuralProbe.__init__(
            self,
            num_channels=num_channels,
            spike_delay=spike_delay,
            spike_peak_duration=spike_peak_duration,
            noise_duration=noise_duration,
            noise_amp_percent=noise_amp_percent,
            fps=fps,
            inner_radius=inner_radius,
            positions_file_path=positions_file_path,
            neighbors_file_path=neighbors_file_path,
            neighbor_radius=neighbor_radius,
            masked_channels=masked_channels
            )

        self.data_file = data_file_path
        positions_file_path = in_probes_dir('positions_mea1k')
        neighbors_file_path = in_probes_dir('neighbormatrix_mea1k')
        if data_file_path is not None:
            d = h5py.File(data_file_path)
            self.d = d
            mapping = d.get('/mapping/')
            channel_indices = mapping['channel']
            electrodes = mapping['electrode'][:]
            routed = np.array(np.where(electrodes > -1))[0]
            self.channels_indices_routed = channel_indices[routed]
            self.nFrames = d['sig'].shape[1] #number_of_frames
            ch_positions = np.vstack((mapping['x'][routed],mapping['y'][routed])).T
            num_channels = ch_positions.shape[0]
            print('# Generating new position and neighbor files from data file')
            create_probe_files(positions_file_path, neighbors_file_path, radius, ch_positions)
        else:
            num_channels = 0
            print('# Note: data file not specified, setting some defaults')

        NeuralProbe.__init__(
            self, num_channels=num_channels, spike_delay=5,
            spike_peak_duration=4, noise_duration=2,
            noise_amp_percent=1, fps=fps,
            inner_radius=100,
            positions_file_path=positions_file_path,
            neighbors_file_path=neighbors_file_path,
            masked_channels=masked_channels)

    def Read(self, t0, t1):
        return self.d['/sig'][self.channels_indices_routed,
                              t0:t1].T.ravel().astype(ctypes.c_short)

class HierlmannVisapyEmulationProbe(NeuralProbe):
    def __init__(self, data_file_path=None, fps=32000, num_channels=102,
                 spike_delay=5, spike_peak_duration=4, noise_duration=3,
                 noise_amp_percent=1, inner_radius=35, neighbor_radius=None,
                 masked_channels=None):
        positions_file_path = in_probe_info_dir('positions_hierlemann_visapy_emulation')
        neighbors_file_path = in_probe_info_dir('neighbormatrix_hierlemann_visapy_emulation')
        NeuralProbe.__init__(
            self,
            num_channels=num_channels,
            spike_delay=spike_delay,
            spike_peak_duration=spike_peak_duration,
            noise_duration=noise_duration,
            noise_amp_percent=noise_amp_percent,
            fps=fps,
            positions_file_path=positions_file_path,
            neighbors_file_path=neighbors_file_path,
            inner_radius=inner_radius, #20
            neighbor_radius=neighbor_radius,
            masked_channels=masked_channels
            )
        self.data_file = data_file_path
        if data_file_path is not None:
            self.d = np.load(self.data_file)
            print('File size: '+str(len(self.d)))
            print('Number of channels: '+str(self.num_channels))
            assert len(self.d) / self.num_channels == len(self.d) // \
                self.num_channels, 'Data not multiple of channel number'
            self.nFrames = len(self.d) // self.num_channels
        else:
            print('# Note: data file not specified, things may break')

    def Read(self, t0, t1):
        return read_flat(self.d, t0, t1, self.num_channels)


class NeuroSeeker_128(NeuralProbe):
    """
     A 128-channel probe designed by the NeuroSeeker project in 2015.
     https://doi.org/10.1109/TRANSDUCERS.2015.7181274
    """
    def __init__(self, data_file_path, num_channels=128, spike_delay=5,
                 spike_peak_duration=5, noise_duration=2, noise_amp_percent=.95,
                 fps=None, inner_radius=1.42, neighbor_radius=None,
                 masked_channels=None):
        positions_file_path = in_probe_info_dir('positions_neuroseeker_128')
        neighbors_file_path = in_probe_info_dir('neighbormatrix_neuroseeker_128')
        NeuralProbe.__init__(
                self,
                num_channels=num_channels,
                spike_delay=spike_delay,
                spike_peak_duration=spike_peak_duration,
                noise_duration=noise_duration,
                noise_amp_percent=noise_amp_percent,
                fps=fps,
                inner_radius=inner_radius,
                positions_file_path=positions_file_path,
                neighbors_file_path=neighbors_file_path,
                neighbor_radius=neighbor_radius,
                masked_channels=masked_channels
                )
        self.data_file = data_file_path
        self.d = openHDF5file(data_file_path)
        self.nFrames, self.fps, self.num_channels, chIndices,\
            self.closest_electrode = getNeuroSeekerParams(self.d, pipette=False)

    def Read(self, t0, t1):
        return readNeuroSeekerProbe(self.d, t0, t1)

class SiNAPS_S1(NeuralProbe):
    """
     A 512-channel probe designed by the NETS3 lab in 2018.
     https://ieeexplore.ieee.org/document/8304606
    """
    def __init__(self, data_file_path, raw_group_path, num_channels=512,
                 spike_delay=5, spike_peak_duration=5, noise_duration=3,
                 noise_amp_percent=.95, fps=25000, inner_radius=40,
                 neighbor_radius=None, masked_channels=None):
        positions_file_path = in_probe_info_dir('positions_SiNAPS_S1')
        neighbors_file_path = in_probe_info_dir('neighbormatrix_SiNAPS_S1')
        NeuralProbe.__init__(
                self,
                num_channels=num_channels,
                spike_delay=spike_delay,
                spike_peak_duration=spike_peak_duration,
                noise_duration=noise_duration,
                noise_amp_percent=noise_amp_percent,
                fps=fps,
                inner_radius=inner_radius,
                positions_file_path=positions_file_path,
                neighbors_file_path=neighbors_file_path,
                neighbor_radius=neighbor_radius,
                masked_channels=masked_channels
                )
        self.data_file = data_file_path
        self.hfile = openHDF5file(data_file_path)
        self.raw_data = self.hfile[raw_group_path]
        self.nFrames = self.raw_data.shape[0]

    def Read(self, t0, t1):
        return readSiNAPS_S1Probe(self.raw_data, t0, t1)
