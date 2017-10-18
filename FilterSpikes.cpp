#include "FilterSpikes.h"

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
			break;
		}
	}
	return areNeighbors;
}

void eliminateDuplicates(Spike max_spike) {

	deque<Spike>::iterator it;
	it = Parameters::spikes_to_be_processed.begin();
	Spike curr_spike;
	int curr_channel, curr_amp;
	
	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		curr_amp = it->amplitude;
		if(areNeighbors(max_spike.channel, curr_channel)) {
			if(max_spike.amplitude > curr_amp) {
				//cout << "Eliminated Dups" << endl;
				//cout << "Spike: " << curr_channel << " " << curr_amp << endl;
				it = Parameters::spikes_to_be_processed.erase(it);
			}
			else {
				++it;
			}
		}
		else {
			++it;
		}
	}
}

Spike filterSpikes(Spike max_spike)
{
	/*Removes all duplicate spikes and returns the max spike

	Parameters
	----------
	max_spike: Spike
		The max candidate for largest amplitude spike.

    Returns
    -------
	max_spike: Spike
		Returns the spike in the neighbors of the max spike that has the
		largest amplitude.
	*/
	int first_spike_channel = max_spike.channel;
	int max_spike_amp = max_spike.amplitude;
	int frame_of_orginal_spike = max_spike.frame;
	int curr_channel, curr_amp, curr_frame;
	//cout << "Original Spike: " << first_spike_channel << " " << max_spike.frame << " " << max_spike_amp << endl;
	Spike curr_spike;
	deque<Spike>::iterator it;
	it = Parameters::spikes_to_be_processed.begin();

	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		curr_amp = it->amplitude;
		curr_frame = it->frame;
		if(areNeighbors(first_spike_channel, curr_channel)) {
			if(curr_frame <= frame_of_orginal_spike + Parameters::noise_duration) {
				if(curr_amp >= max_spike_amp) {
					it = Parameters::spikes_to_be_processed.erase(it);
					max_spike = curr_spike;
					max_spike_amp = curr_amp;
				}
				else {
					it = Parameters::spikes_to_be_processed.erase(it);
				}
			}
			else {
				if(curr_amp >= max_spike_amp*Parameters::noise_amp_percent) {
					++it;
				}
				else {
					it = Parameters::spikes_to_be_processed.erase(it);
				}
			}
		}
		else {
			++it;
		}
	}

	if(Parameters::spikes_to_be_processed.size() != 0) {
		eliminateDuplicates(max_spike);
	}

	return max_spike;
}

}
