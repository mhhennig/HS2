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

using namespace std;

//Internal representation of a spike. User has no need to use it.
struct Spike {
	int amplitude;
	int channel;
	int frame;
	vector<int> amp_cutouts;
	vector<int> written_cutout;
};

namespace Parameters {

extern int num_channels; //Number of channels on the probe
extern int num_recording_channels; //Number of channels to be used for spike data
extern int spike_delay; //The number of frames back a spike occurred after it was detected (where the beginning of the spike was).
extern int spike_peak_duration; //The number of frames it takes a spike amplitude to fully decay.
extern int noise_duration; //The number of frames that the true spike can occur after the first detection.
extern int noise_amp; //Amplitude difference allowed to differentiate between increasing amplitude duplicates
extern int max_neighbors;//Maximum number of neighbors a channel can have in the probe
extern int** neighbor_matrix;/*Indexed by the channel number starting at 0 and going up to num_recording_channels - 1. Each 
							  index contains pointer to another array which contains channel number of all its neighbors.
							  User creates this before calling SpikeHandler. Each column has size equal to max neighbors where
							  any channels that have less neighbors fills the rest with -1 (important). */
extern int** channel_positions;/*Indexed by the channel number starting at 0 and going up to num_recording_channels - 1. Each 
							  index contains pointer to another array which contains X and Y position of the channel. User creates
							  this before calling SpikeHandler. */
extern int aGlobal; //Global noise
extern int** baselines; //Contains spike_delay number of frames of median baseline values. Updated by user at every frame.
extern bool to_localize; //True: filter and localize the spike, False: just filter the spike.
extern deque<Spike> spikes_to_be_processed; //Contains all spikes to be proccessed when spike_peak_duration number of frames is stored.
extern int cutout_length; //The cutout length to be written out for the spike. Can't be larger than extra data tacked on to raw data.
extern int filtered_spikes; //number of filtered spikes
extern short* raw_data; //raw data passed in for current iteration
extern int index_data; //The index given to start accessing the raw data. To account for extra data tacked on for cutout purposes.
extern int index_baselines; /*The index given to start accessing the baseline array since baseline array is size 5 and location of
							  oldest baseline is constantly changing*/
extern int frames; //Number of current iterations of raw data passed in. User starts this at 0 and increments it for each chunk of data;
extern int iterations; //The number of frames passed into loadRawData EXCLUDING the buffer frames.
extern int maxsl; //Number of frames after a detection that a spike is accepted

};

#endif