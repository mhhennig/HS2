from __future__ import division
import numpy as np

# import json
from matplotlib import pyplot as plt
from .probe_functions.readUtils import read_flat
from .probe_functions.readUtils import openHDF5file, getHDF5params, getHDF5params_brw4
from .probe_functions.readUtils import readHDF5t_100, readHDF5t_101
from .probe_functions.readUtils import readHDF5t_100_i, readHDF5t_101_i
from .probe_functions.readUtils import readHDF5_brw4
from .probe_functions.neighborMatrixUtils import createNeighborMatrix

# import h5py
import ctypes
import os.path
from scipy.spatial.distance import cdist
import warnings

this_file_path = "."  # os.path.dirname(os.path.abspath(__file__))

DEFAULT_EVENT_LENGTH = 0.5
DEFAULT_PEAK_JITTER = 0.2


def get_neighbors(radius, ch_positions):
    distances = cdist(ch_positions, ch_positions, metric="euclidean")
    indices = np.arange(ch_positions.shape[0])
    neighbors = []
    for dist_from_ch in distances:
        neighbors.append(indices[dist_from_ch <= radius])
    return neighbors


# def create_probe_files(pos_file, neighbor_file, radius, ch_positions):
#     n_channels = ch_positions.shape[0]
#     # NB: Notice the column, row order in write
#     with open(pos_file, "w") as f:
#         for pos in ch_positions:
#             f.write("{},{},\n".format(pos[0], pos[1]))
#     f.close()
#     # # NB: it is also possible to use metric='cityblock' (Manhattan distance)
#     distances = cdist(ch_positions, ch_positions, metric="euclidean")
#     indices = np.arange(n_channels)
#     with open(neighbor_file, "w") as f:
#         for dist_from_ch in distances:
#             neighbors = indices[dist_from_ch <= radius]
#             f.write("{},\n".format(str(list(neighbors))[1:-1]))
#     f.close()


def in_probes_dir(file):
    probe_path1 = os.getenv("HS2_PROBE_PATH", this_file_path)
    probe_path = os.path.join(probe_path1, "probes")
    if not os.path.exists(probe_path):
        os.mkdir(probe_path)
    return os.path.join(probe_path, file)


def in_probe_info_dir(file):
    probe_path1 = os.getenv("HS2_PROBE_PATH", this_file_path)
    probe_path = os.path.join(probe_path1, "probe_info")
    if not os.path.exists(probe_path):
        os.mkdir(probe_path)
    return os.path.join(probe_path, file)


class NeuralProbe(object):
    def __init__(
        self,
        num_channels,
        noise_amp_percent,
        inner_radius,
        fps,
        positions_file_path,
        neighbors_file_path,
        neighbor_radius,
        event_length,
        peak_jitter,
        masked_channels=[],
        spike_peak_duration=None,
        noise_duration=None,
    ):
        if neighbor_radius is not None and positions_file_path is not None:
            createNeighborMatrix(
                neighbors_file_path, positions_file_path, neighbor_radius
            )
        self.fps = fps
        self.num_channels = num_channels
        self.spike_peak_duration = self._deprecate_or_convert(
            spike_peak_duration, event_length, "spike_peak_duration", "event_length"
        )
        self.noise_duration = self._deprecate_or_convert(
            noise_duration, peak_jitter, "noise_duration", "peak_jitter"
        )
        self.noise_amp_percent = noise_amp_percent
        self.positions_file_path = positions_file_path
        self.neighbors_file_path = neighbors_file_path
        self.masked_channels = masked_channels
        self.inner_radius = inner_radius
        if masked_channels is None:
            self.masked_channels = []

        if positions_file_path is not None and neighbors_file_path is not None:
            self.loadPositions(positions_file_path)
            self.loadNeighbors(neighbors_file_path)

    # Load in neighbor and positions files
    def loadNeighbors(self, neighbors_file_path):
        neighbor_file = open(neighbors_file_path, "r")
        neighbors = []
        for neighbor in neighbor_file.readlines():
            neighbors.append(np.array(neighbor[:-2].split(",")).astype(int))
        neighbor_file.close()
        # assert len(neighbors) == len(pos)
        self.neighbors = neighbors
        self.max_neighbors = max([len(n) for n in neighbors])

    def loadPositions(self, positions_file_path):
        position_file = open(positions_file_path, "r")
        positions = []
        for position in position_file.readlines():
            positions.append(np.array(position[:-2].split(",")).astype(float))
        self.positions = np.asarray(positions)
        position_file.close()

    def _deprecate_or_convert(self, old_var, new_var, old_name, new_name):
        if old_var is not None:
            warnings.warn(
                "{} is deprecated and will be removed. ".format(old_name)
                + "Set {} instead (in milliseconds). ".format(new_name)
                + "{} takes priority over {}!".format(old_name, new_name),
                DeprecationWarning,
            )
            return int(old_var)
        else:
            return int(new_var * self.fps / 1000)

    # Show visualization of probe
    def show(self, show_neighbors=[10], figwidth=3):
        xmax, ymax = self.positions.max(0)
        xmin, ymin = self.positions.min(0)
        ratio = ymax / xmax
        plt.figure(figsize=(figwidth, figwidth * ratio))
        for ch in show_neighbors:
            for neighbor in self.neighbors[ch]:
                plt.plot(
                    [self.positions[ch, 0], self.positions[neighbor, 0]],
                    [self.positions[ch, 1], self.positions[neighbor, 1]],
                    "--k",
                    alpha=0.7,
                )
        plt.scatter(*self.positions.T)
        plt.scatter(*self.positions[self.masked_channels].T, c="r")
        for i, pos in enumerate(self.positions):
            plt.annotate(i, pos)

    def Read(self, t0, t1):
        raise NotImplementedError(
            "The Read function is not implemented for \
            this probe"
        )

    def getChannelsPositions(self, channels):
        channel_positions = []
        for channel in channels:
            if channel >= self.num_channels:
                raise ValueError(
                    "Channel index too large, maximum " + self.num_channels
                )
            else:
                channel_positions.append(self.positions[channel])
        return channel_positions


class BioCam(NeuralProbe):
    def __init__(
        self,
        data_file_path=None,
        num_channels=4096,
        fps=0,
        noise_amp_percent=1,
        inner_radius=1.75,
        neighbor_radius=None,
        masked_channels=[0],
        event_length=DEFAULT_EVENT_LENGTH,
        peak_jitter=DEFAULT_PEAK_JITTER,
    ):
        self.data_file = data_file_path
        if data_file_path is not None:
            self.d = openHDF5file(data_file_path)
            if "3BData" in self.d:
                params = getHDF5params(self.d)
            else:
                params = getHDF5params_brw4(self.d)
            self.nFrames, sfd, self.num_channels, chIndices, file_format, inversion = (
                params
            )
            print(
                "# Signal inversion looks like",
                inversion,
                ", guessing the "
                "right method for data access.\n# If your detection results "
                "look strange, signal polarity is wrong.\n",
            )
            if file_format == 100:
                if inversion == -1:
                    self.read_function = readHDF5t_100
                else:
                    self.read_function = readHDF5t_100_i
            elif file_format == 101 or file_format == 102:
                if inversion == -1:
                    self.read_function = readHDF5t_101_i
                else:
                    self.read_function = readHDF5t_101
            elif file_format == "brw4":
                self.read_function = readHDF5_brw4
        else:
            print("# Note: data file not specified, setting some defaults")
            self.num_channels = 4096
            sfd = fps
        if self.num_channels < 4096:
            print(
                "# Note: only",
                self.num_channels,
                "channels recorded, fixing positions/neighbors",
            )
            print("# This may break - known to work only for rectangular sections!")
            recorded_channels = self.d["3BRecInfo"]["3BMeaStreams"]["Raw"]["Chs"]
        else:
            recorded_channels = None

        positions_file_path = in_probe_info_dir("positions_biocam")
        neighbors_file_path = in_probe_info_dir("neighbormatrix_biocam")

        NeuralProbe.__init__(
            self,
            num_channels=self.num_channels,
            noise_amp_percent=noise_amp_percent,
            fps=sfd,
            inner_radius=inner_radius,  # 1.75,
            positions_file_path=positions_file_path,
            neighbors_file_path=neighbors_file_path,
            neighbor_radius=neighbor_radius,
            masked_channels=masked_channels,
            event_length=event_length,
            peak_jitter=peak_jitter,
        )

        # if a probe only records a subset of channels, filter out unused ones
        # this may happen in Biocam recordings
        # note this uses positions to identify recorded channels
        # requires channels are an ordered list
        if recorded_channels is not None:
            inds = np.zeros(self.num_channels, dtype=int)
            for i, c in enumerate(recorded_channels):
                inds[i] = np.where(
                    np.all(
                        (self.positions - np.array([c[1] - 1, c[0] - 1])) == 0, axis=1
                    )
                )[0]
            self.positions = self.positions[inds].astype(int)
            x0 = np.min([p[0] for p in self.positions])
            y0 = np.min([p[1] for p in self.positions])
            x1 = np.max([p[0] for p in self.positions])
            y1 = np.max([p[1] for p in self.positions])
            print("# Array boundaries (x):", x0, x1)
            print("# Array boundaries (y):", y0, y1)
            print("# Array width and height:", x1 - x0 + 1, y1 - y0 + 1)
            print("# Number of channels:", self.num_channels)
            lm = np.zeros((64, 64), dtype=int) - 1
            # oddness because x/y are transposed in brw
            lm[y0 : y1 + 1, x0 : x1 + 1] = np.arange(self.num_channels).T.reshape(
                y1 - y0 + 1, x1 - x0 + 1
            )
            self.neighbors = [lm.flatten()[self.neighbors[i]] for i in inds]
            self.neighbors = [n[(n >= 0)] for n in self.neighbors]
            self.positions = self.positions - np.min(self.positions, axis=0)

    def Read(self, t0, t1):
        return self.read_function(self.d, t0, t1, self.num_channels)


class RecordingExtractor(NeuralProbe):
    def __init__(
        self,
        re,
        noise_amp_percent=1,
        inner_radius=60,
        neighbor_radius=60,
        masked_channels=None,
        xy=None,
        event_length=DEFAULT_EVENT_LENGTH,
        peak_jitter=DEFAULT_PEAK_JITTER,
    ):
        self.d = re
        # positions_file_path = in_probe_info_dir("positions_spikeextractor")
        # neighbors_file_path = in_probe_info_dir("neighbormatrix_spikeextractor")
        try:
            self.nFrames = re.get_num_frames()
        except:
            self.nFrames = re.get_num_frames(0)
        num_channels = re.get_num_channels()
        fps = re.get_sampling_frequency()
        ch_positions = np.array(
            [
                np.array(re.get_channel_property(ch, "location"))
                for ch in re.get_channel_ids()
            ]
        )
        if ch_positions.shape[1] > 2:
            if xy is None:
                print(
                    "# Warning: channel locations have",
                    ch_positions.shape[1],
                    "dimensions",
                )
                print("#          using the last two.")
                xy = (ch_positions.shape[1] - 2, ch_positions.shape[1] - 1)
            ch_positions = ch_positions[:, xy]
        print("# Generating new position and neighbor files from data file")
        # create_probe_files(
        #     positions_file_path, neighbors_file_path, inner_radius, ch_positions
        # )
        self.positions = ch_positions
        self.neighbors = get_neighbors(inner_radius, ch_positions)
        self.max_neighbors = max([len(n) for n in self.neighbors])

        NeuralProbe.__init__(
            self,
            num_channels=num_channels,
            noise_amp_percent=noise_amp_percent,
            fps=fps,
            inner_radius=inner_radius,
            positions_file_path=None,
            neighbors_file_path=None,
            masked_channels=masked_channels,
            neighbor_radius=neighbor_radius,
            event_length=event_length,
            peak_jitter=peak_jitter,
        )

    def Read(self, t0, t1):
        return (
            self.d.get_traces(
                channel_ids=self.d.get_channel_ids(), start_frame=t0, end_frame=t1
            )
            .ravel()
            .astype(ctypes.c_short)
        )
