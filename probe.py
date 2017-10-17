import numpy as np
from matplotlib import pyplot as plt


class NeuralProbe(object):
    def __init__(self, num_channels, spike_delay,
                 spike_peak_duration, noise_duration, noise_amp,
                 fps, positions_path, neighbors_path):
        self.num_channels = num_channels
        self.spike_delay = spike_delay
        self.spike_peak_duration = spike_peak_duration
        self.noise_duration = noise_duration
        self.noise_amp = noise_amp
        self.fps = fps
        self.positions_path = positions_path
        self.neighbors_path = neighbors_path
        f = open(neighbors_path, 'r')
        neighs = []
        for l in f.readlines():
            neighs.append(np.array(l[:-2].split(',')).astype(int))
        f.close()
        f = open(positions_path, 'r')
        pos = []
        for l in f.readlines():
            pos.append(np.array(l[:-2].split(',')).astype(int))
        f.close()
        assert len(neighs) == len(pos)
        self.num_recording_channels = len(neighs)
        self.positions = np.asarray(pos)
        self.neighbors = neighs
        self.max_neighbors = max([len(n) for n in neighs])

    def Show(self, showNeighbors=[10], figwidth=3):
        xmax, ymax = self.positions.max(0)
        xmin, ymin = self.positions.min(0)
        ratio = ymax/xmax
        plt.figure(figsize=(figwidth, figwidth*ratio))
        for ch in showNeighbors:
            for neigh in self.neighbors[ch]:
                plt.plot([self.positions[ch, 0], self.positions[neigh, 0]],
                         [self.positions[ch, 1], self.positions[neigh, 1]],
                         '--k', alpha=.7)
        plt.scatter(*self.positions.T)
        for i, pos in enumerate(self.positions):
            plt.annotate(i, pos)
        plt.ylim([0, ymax+ymin])
        plt.xlim([0, xmax+xmin])
