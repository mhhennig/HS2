#ifndef PARAMETERS_H
#define PARAMETERS_H

//Contains all parameters and libraries for running the SpikeHandler Methods

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
#include <iterator>
#include <vector>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>
//#include<bits/stdc++.h>

using namespace std;

//Internal representation of a spike. User has no need to use it.
struct Spike {
	int amplitude;
	int channel;
	int frame;
	deque<int> largest_channels;
	vector<int32_t> written_cutout;
	tuple<vector<int>,int*> waveformscounts;
    // //These contain all information of what occurred at neighbors
    // vector<Event> inner_neighbors;
    // vector<Event> outer_neighbors;
};

namespace Parameters {

extern int ASCALE; //Scaling on the raw extracellular data
extern int num_com_centers; //Number of channels used for center of mass
extern int num_channels; //Number of channels on the probe
extern int spike_peak_duration; //The number of frames it takes a spike amplitude to fully decay.
extern int noise_duration; //The number of frames that the true spike can occur after the first detection.
extern float noise_amp_percent; //Amplitude percentage allowed to differentiate between decreasing amplitude duplicate spike
extern int max_neighbors;//Maximum number of neighbors a channel can have in the probe
extern int** neighbor_matrix;/*Indexed by the channel number starting at 0 and going up to num_channels - 1. Each
							  index contains pointer to another array which contains channel number of all its neighbors.
							  User creates this before calling SpikeHandler. Each column has size equal to max neighbors where
							  any channels that have less neighbors fills the rest with -1 (important). */
extern int** inner_neighbor_matrix; /*Indexed by the channel number starting at 0 and going up to num_channels - 1. Each
							  index contains pointer to another array which contains channel number of all its inner neighbors.
							  Created by SpikeHandler; */
extern int** outer_neighbor_matrix; /*Indexed by the channel number starting at 0 and going up to num_channels - 1. Each
                                    index contains pointer to another array which contains channel number of all its outer neighbors.
                              		Created by SpikeHandler; */

extern float** channel_positions;/*Indexed by the channel number starting at 0 and going up to num_channels - 1. Each
							  index contains pointer to another array which contains X and Y position of the channel. User creates
							  this before calling SpikeHandler. */
extern int aGlobal; //Global noise
extern int** baselines; //Contains spike_delay number of frames of median baseline values. Updated by user at every frame.
extern bool to_localize; //True: filter and localize the spike, False: just filter the spike.
extern deque<Spike> spikes_to_be_processed; //Contains all spikes to be proccessed when spike_peak_duration number of frames is stored.
extern int cutout_start; //The number of frames before the spike that the cutout starts at
extern int cutout_end; //The number of frames after the spike that the cutout ends atextern int filtered_spikes; //number of filtered spikes
extern short* raw_data; //raw data passed in for current iteration
extern int index_data; //The index given to start accessing the raw data. To account for extra data tacked on for cutout purposes.
extern int index_baselines; /*The index given to start accessing the baseline array since baseline array is size 5 and location of
							  oldest baseline is constantly changing*/
extern int max_frames_processed; //The number of frames passed into loadRawData EXCLUDING the buffer frames.
extern int before_chunk; //The number of buffer frames before chunk
extern int after_chunk; //The number of buffer frames after chunk
extern int iterations; //Number of current iterations of raw data passed in. User starts this at 0 and increments it for each chunk of data;
extern int maxsl; //Number of frames after a detection that a spike is accepted
extern int end_raw_data; //index of the end of the raw data
extern int* masked_channels; //stores all masked channels as 0 and regular channels as 1
extern int event_number;
extern float inner_radius;
extern bool debug;
extern bool verbose;
extern bool decay_filtering; //if true, then tries to filter by decay (more effective for less dense arrays)

};

#endif
