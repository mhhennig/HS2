#include "FilterSpikes.h"
#include <iostream>

using namespace std;

namespace FilterSpikes {
	
bool areNeighbors(int channel_one, int channel_two)
{
	/*Determines whether or not two channels are neighbors.

	Parameters
	----------
	channel_one: int
		The channel number whose neighborhoold is checked.
	channel_two: int
		The channel number of the potential neighbor.

	Returns
	-------
	isNeighbor: bool
		True if the channels are neighbors.
		False if the channels are not neighbors.
	*/
	bool areNeighbors = false;
	for(int i = 0; i < Parameters::max_neighbors; i++) {
		if(Parameters::neighbor_matrix[channel_one][i] == channel_two) {
			areNeighbors = true;
		}
	}
	return areNeighbors;
}

Spike filterSpikes(Spike curr_original_spike)
{
	/*Removes all duplicate spikes and returns the original spike

	Parameters
	----------
	curr_original_spike: Spike
		The original candidate for largest amplitude spike.
    Returns
    -------
	original_spike: Spike
		Returns the spike in the neighbors of the original spike that has the
		largest amplitude.
	*/
	int first_spike_channel = curr_original_spike.channel;
	int curr_original_spike_amp = curr_original_spike.amplitude;
	int frame_of_orginal_spike = curr_original_spike.frame;
	deque<Spike>::iterator it;

	it = Parameters::spikes_to_be_processed.begin();
	int curr_channel, curr_amp, curr_frame;
	Spike curr_spike;
	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		curr_amp = it->amplitude;
		curr_frame = it->frame;
		//cout << "Are Neighbors" << '\n';
		bool channelsAreNeighbors = areNeighbors(first_spike_channel, curr_channel);
		//cout << "Are neighbors OK" << '\n';
		if(channelsAreNeighbors && curr_amp >= curr_original_spike_amp + Parameters::noise_amp) {
			if(curr_frame <= frame_of_orginal_spike + Parameters::noise_duration) {
				it = Parameters::spikes_to_be_processed.erase(it);
				curr_original_spike = curr_spike;
				curr_original_spike_amp = curr_amp;
			}
			else {
				++it;
			}
		}
		else if(channelsAreNeighbors) {
			it = Parameters::spikes_to_be_processed.erase(it);
		}
		else {
			++it;
		}
	}
	return curr_original_spike;
}
}