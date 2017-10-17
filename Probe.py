import numpy as np
from matplotlib import pyplot as plt


class NeuralProbe(object):
    def __init__(self, num_channels, spike_delay,
                 spike_peak_duration, noise_duration, noise_amp,
                 fps):
        self.num_channels = num_channels
        self.spike_delay = spike_delay
        self.spike_peak_duration = spike_peak_duration
        self.noise_duration = noise_duration
        self.noise_amp = noise_amp
        self.fps = fps
        self.num_recording_channels = None;
        self.positions = None;
        self.neighbors = None;
        self.max_neighbors = None;

    #Parameters they may want to change after testing on the data
    def setSpikeDelay(self, spike_delay):
    	self.spike_delay = spike_delay

    def setSpikePeakDuration(self, spike_peak_duration):
    	self.spike_peak_duration = spike_peak_duration

    def setNoiseDuration(self, noise_duration):
    	self.noise_duration = noise_duration

    def setNoiseAmp(self, noise_amp):
    	self.noise_amp = noise_amp

    #Parameters they may be interested in looking at outside of the detection	
    def getSpikeDelay(self):
    	return self.spike_delay

    def getSpikePeakDuration(self):
    	return self.spike_peak_duration

    def getNoiseDuration(self):
    	return self.noise_duration

    def getNoiseAmp(self):
    	return self.noise_amp

    def getFps(self):
    	return self.fps

    def getNumChannels(self):
        if self.num_channels == None:
            print("num_channels not set yet")
        return self.num_channels

    def getNumRecordingChannels(self):
    	if self.num_recording_channels == None:
    		print("num_recording_channels not set yet")
    	return self.num_recording_channels

    def getPositions(self):
    	if type(self.positions) is not np.ndarray:
    		print("positions not set yet")
    	return self.positions

    def getNeighbors(self):
    	if type(self.neighbors) is not np.ndarray:
    		print("neighbors not set yet")
    	return self.neighbors

    def getMaxNeighbors(self):
    	if self.max_neighbors == None:
    		print("max_neighbors not set yet")
    	return self.max_neighbors

    #Load in neighbor and positions files
    def loadNeighbors(self, neighbors_file_path):
        neighbor_file = open(neighbors_file_path, 'r')
        neighbors = []
        for neighbor in neighbor_file.readlines():
            neighbors.append(np.array(neighbor[:-2].split(',')).astype(int))
        neighbor_file.close()
        #assert len(neighbors) == len(pos)
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

    #Show visualization of probe
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

    def __init__(self):
        NeuralProbe.__init__(self, num_channels=385, spike_delay=5,
        					 spike_peak_duration=5, noise_duration=3,
        					 noise_amp = 80000, fps=30000)
       	self.loadNeighbors("neighbormatrix")
       	self.loadPositions('positions')
