#ifndef SPIKEHANDLER_H  
#define SPIKEHANDLER_H 

#include "ProcessSpikes.h"

using namespace std;

void setInitialParameters(int _num_channels, int _num_recording_channels, int _spike_delay, int _spike_peak_duration, \
						  int _noise_duration, int _noise_amp, int** _channel_positions, int** _neighbor_matrix, \
						  int _max_neighbors, bool _to_localize, int _cutout_length);
void loadRawData(short* _raw_data, int iterations, int frames, int _index_data);
void setLocalizationParameters(int _aGlobal, int** _baselines);
void addSpike(int channel, int frame, int amplitude);
void terminateSpikeHandler();


#endif