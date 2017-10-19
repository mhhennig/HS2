import numpy as np
from matplotlib import pyplot as plt


class NeuralProbe(object):
    def __init__(self, num_channels, spike_delay,
                 spike_peak_duration, noise_duration, noise_amp_percent,
                 fps, positions_file_path, neighbors_file_path):
        self.num_channels = num_channels
        self.spike_delay = spike_delay
        self.spike_peak_duration = spike_peak_duration
        self.noise_duration = noise_duration
        self.noise_amp_percent = noise_amp_percent
        self.fps = fps
        self.positions_file_path = positions_file_path
        self.neighbors_file_path = neighbors_file_path

        self.loadPositions(positions_file_path)
        self.loadNeighbors(neighbors_file_path)

    # Load in neighbor and positions files
    def loadNeighbors(self, neighbors_file_path):
        neighbor_file = open(neighbors_file_path, 'r')
        neighbors = []
        for neighbor in neighbor_file.readlines():
            neighbors.append(np.array(neighbor[:-2].split(',')).astype(int))
        neighbor_file.close()
        # assert len(neighbors) == len(pos)
        self.num_recording_channels = len(neighbors)
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
        ratio = ymax/xmax
        plt.figure(figsize=(figwidth, figwidth*ratio))
        for ch in show_neighbors:
            for neighbor in self.neighbors[ch]:
                plt.plot([self.positions[ch, 0], self.positions[neighbor, 0]],
                         [self.positions[ch, 1], self.positions[neighbor, 1]],
                         '--k', alpha=.7)
        plt.scatter(*self.positions.T)
        for i, pos in enumerate(self.positions):
            plt.annotate(i, pos)
        plt.ylim([0, ymax+ymin])
        plt.xlim([0, xmax+xmin])


class NeuroPixel(NeuralProbe):
    def __init__(self, fps=30000):
        NeuralProbe.__init__(self, num_channels=385, spike_delay=5,
                             spike_peak_duration=5, noise_duration=2,
                             noise_amp_percent=.95, fps=fps,
                             positions_file_path='positions',
                             neighbors_file_path='neighbormatrix')

class BioCam(NeuralProbe):
    def __init__(self, fps=0):
        NeuralProbe.__init__(self, num_channels=4096, spike_delay=5,
                             spike_peak_duration=5, noise_duration=2,
                             noise_amp_percent=.95, fps=fps,
                             positions_file_path='positions_biocam',
                             neighbors_file_path='neighbormatrix_biocam')
