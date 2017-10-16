import numpy as np
import pandas as pd
from detect import detectData


class herdingspikes(object):
    def __init__(self, classpath=None, datapath=None):
        if classpath is not None:
            #  we want this to load a preexisting instantiation
            #  of detected stuff
            raise NotImplementedError()

    def DetectFromRaw(self, num_channels, num_recording_channels, spike_delay,
                      spike_peak_duration, noise_duration, noise_amp,
                      max_neighbors, to_localize, cutout_length, fps, threshold,
                      maa=0, maxsl=12, minsl=3, ahpthr=0):
        detectData(self.datapath, num_channels, num_recording_channels,
                   spike_delay, spike_peak_duration, noise_duration, noise_amp,
                   max_neighbors, to_localize, cutout_length, fps, threshold,
                   maa, maxsl, minsl, ahpthr)
        sp = np.loadtxt('ProcessedSpikes')
        self.spikes = pd.DataFrame({'Channel': sp[:, 0].astype(int)
                                    'Time': sp[:, 1].astype(int)
                                    'Amplitude': sp[:, 2]
                                    'x': sp[:, 3]
                                    'y': sp[:, 4]
                                    })
