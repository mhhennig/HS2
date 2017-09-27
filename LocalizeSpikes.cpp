#include "LocalizeSpikes.h"

tuple<float, float> tuple<float,float> localizeSpike(Spike spike_to_be_localized)
{
	/*Estimates the X and Y position of where a spike occured on the probe.

	Parameters
	----------
	neighbor_matrix: np.array
		A array with each index representing channel
		numbers that correspond to integer array values that contain the channel
		numbers that the index channel is neighboring.
	spike_to_be_localized: list
		The spike that will be used to determine where the spike occurred. Each
		spike is represented by a list containing channel number, frame, and amplitude.
	baselines: dict
		A dictionary of dictionaries that contains all frames up to 1000
		and at each frame, every channel stores a baseline value
	aGlobals: dict
		A dictionary of frames with their corresponding noise amplitudes
	raw_data: numpy.core.memmap.memmap
		Contains all raw data from the sensors
	channel_positions: list
		A list of tuples that stores the X and Y positions of every channel
	wave_delay:
		How many frames back the spike occurs after the frame the spike is detected.
	Returns
	-------
	position: tuple
		An X and Y coordinate tuple that corresponds to where the spike occurred
	*/
	int ASCALE = -64;
	int spike_channel = spike_to_be_localized.channel;
	int spike_frame = spike_to_be_localized.frame;
	int spike_start_frame = spike_frame - wave_delay;
	int spike_end_frame = spike_frame + wave_delay;
	tuple<int,int> channel_amps[Parameters::max_neighbors];
	int curr_largest_amp = 10000;
	int curr_reading, curr_amp, start_cutout, end_cutout, curr_neighbor_channel;

	for(int i = 0; i < MAX_NEIGHBORS; i++) {
		curr_neighbor_channel = neighbor_matrix[spike_channel][i];
		//Out of neighbors
		if(curr_neighbor_channel = -1) {
			break;
		}
		start_cutout = spike_start_frame*MAX_CHANNELS + curr_neighbor_channel;
		end_cutout = spike_end_frame*MAX_CHANNELS + curr_neighbor_channel;
		for(int j = 0; j < end_cutout; j++) {
			curr_reading = raw_data[start_cutout + j];
			curr_amp = curr_reading - aGlobal * ASCALE - baselines[curr_neighbor_channel];
			if(curr_amp < curr_largest_amp) {
				curr_largest_amp = curr_amp;
			}
		}
		channel_amps[i] = make_tuple(curr_neighbor_channel, curr_largest_amp);
	}

	tuple<float,float> centerOfMass = centerOfMass(channel_amps, channel_positions);
	return centerOfMass;
}