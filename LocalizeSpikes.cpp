#include "LocalizeSpikes.h"

namespace LocalizeSpikes {

//tuple<int,int> centerOfMass(int spike_channel)
int centerOfMass(int spike_channel)
{
	int X = 0;
	int Y = 0;
	int X_numerator = 0;
	int Y_numerator = 0;
	int denominator = 0;
	int X_coordinate;
	int Y_coordinate;
	int channel;
	int weight;
	int amps_size = Parameters::amps.size();

	for(int i = 0; i < amps_size; i++) {
		//cout << amps->at(i) << '\n';
		if(Parameters::amps.at(i) > 0) {
			channel = Parameters::neighbor_matrix[spike_channel][i];
			X_coordinate = Parameters::channel_positions[channel][0];
			Y_coordinate = Parameters::channel_positions[channel][1];
			weight = Parameters::amps.at(i);
			X_numerator += weight * X_coordinate;
			Y_numerator += weight * Y_coordinate;
			denominator += weight;
		}
	}

	X = X_numerator / denominator;
	Y = Y_numerator / denominator;

	//return make_tuple(X, Y);
	return X;
}

//tuple<int, int> localizeSpike(Spike spike_to_be_localized)
int localizeSpike(Spike spike_to_be_localized)
{
	/*Estimates the X and Y position of where a spike occured on the probe.

	Parameters
	----------
	spike_to_be_localized: Spike
		The spike that will be used to determine where the origin of the spike occurred.
	Returns
	-------
	position: tuple<float, float>
		An X and Y coordinate tuple that corresponds to where the spike occurred.
	*/
	int ASCALE = -64;
	int spike_channel = spike_to_be_localized.channel;
	int spike_frame = spike_to_be_localized.frame;
	int spike_start_frame = spike_frame - Parameters::spike_delay;
	int spike_end_frame = spike_frame + Parameters::spike_delay;
	int curr_largest_amp = -100000; //arbitrarily small to make sure that it is immediately overwritten
	int curr_reading, start_cutout, end_cutout, curr_neighbor_channel;
	int curr_amp;

	for(int i = 0; i < Parameters::max_neighbors; i++) {
		curr_neighbor_channel = Parameters::neighbor_matrix[spike_channel][i];
		if(curr_neighbor_channel != -1) {
			start_cutout = spike_start_frame*Parameters::num_channels + curr_neighbor_channel;
			end_cutout = spike_end_frame*Parameters::num_channels + curr_neighbor_channel;
			for(int j = 0; j < end_cutout; j++) {
				curr_reading = Parameters::raw_data[start_cutout + j];
				curr_amp = -1*(curr_reading - Parameters::aGlobal * ASCALE - Parameters::baselines[curr_neighbor_channel]);
				if(curr_amp > curr_largest_amp) {
					curr_largest_amp = curr_amp;
				}
			}
			Parameters::amps.push_back(curr_largest_amp);
			curr_largest_amp = -10000;
		}
	}

	sort(begin(Parameters::amps), end(Parameters::amps)); //sort the array
	
	//Find median of array
	int amps_size = Parameters::amps.size();
	int median;
	if(amps_size % 2 == 0) {
		median = (Parameters::amps.at(amps_size/2) + Parameters::amps.at(amps_size/2 + 1))/2;
	}
	else {
		median = Parameters::amps.at(amps_size/2);	
	}

	//Center amplitudes
	for(int i = 0; i < amps_size; i++) {
		Parameters::amps.at(i) = Parameters::amps.at(i) - median;
	}

	//tuple<int,int> position = centerOfMass(spike_channel);
	int position = centerOfMass(spike_channel);
	Parameters::amps.clear();

	return position;
}
}