#ifndef PARAMETERS_H  
#define PARAMETERS_H 

#include <deque>
#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <deque>
#include <tuple>

using namespace std;

struct Spike {
	int amplitude;
	int channel;
	int frame;
};

namespace Parameters {

extern int num_channels;
extern int num_recording_channels;
extern int spike_delay;
extern int spike_peak_duration;
extern int noise_duration;
extern int noise_amp;
extern int max_neighbors;
extern int** neighbor_matrix;
extern int** channel_positions;
extern int aGlobal;
extern int* baselines;
extern bool to_localize;
extern deque<Spike> spikes_to_be_processed;
extern int start_cutout;
extern int end_cutout;
extern int filtered_spikes;
extern short* raw_data;
extern deque<int> amps;

};

#endif