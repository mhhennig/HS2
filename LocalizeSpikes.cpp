#include "LocalizeSpikes.h"

namespace LocalizeSpikes {

tuple<int, int> centerOfMass(int spike_channel)
{
	/*Calculates the center of mass of a spike to calculate where it occurred using a weighted average.

	Parameters
	----------
	spike_channel: int
		The channel at which the true spike was detected.

	Returns
	-------
	position: tuple<float, float>
		An X and Y coordinate tuple that corresponds to where the spike occurred.
	*/
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

	return make_tuple(X, Y);
}

tuple<int, int> localizeSpike(Spike spike_to_be_localized, int baseline_frame)
{
	/*Estimates the X and Y position of where a spike occured on the probe.

	Parameters
	----------
	spike_to_be_localized: Spike
		The spike that will be used to determine where the origin of the spike occurred.
	baseline_frame: int
		The frame difference from the first detected spike that the true spike is detected.
		Used to access the baseline array.

	Returns
	-------
	position: tuple<float, float>
		An X and Y coordinate tuple that corresponds to where the spike occurred.
	*/
	int ASCALE = -64;
	int spike_channel = spike_to_be_localized.channel;
	int curr_largest_amp = -100000; //arbitrarily small to make sure that it is immediately overwritten
	int curr_reading, curr_neighbor_channel;
	int curr_amp;

	int amp_cutout_size = spike_to_be_localized.amp_cutouts.size();
	for(int i = 0; i < amp_cutout_size; i++) {
		curr_neighbor_channel = Parameters::neighbor_matrix[spike_channel][i/(Parameters::spike_delay*2)];
		curr_reading = spike_to_be_localized.amp_cutouts.at(i);
		curr_amp = -1*(curr_reading - Parameters::aGlobal * ASCALE - Parameters::baselines[curr_neighbor_channel][(Parameters::index_baselines - baseline_frame) % (Parameters::spike_delay + 1)]);
		if(curr_amp > curr_largest_amp) {
			curr_largest_amp = curr_amp;
		}
		if(i % (Parameters::spike_delay*2) == Parameters::spike_delay*2 - 1) {
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
	
	tuple<int,int> position = centerOfMass(spike_channel);
	Parameters::amps.clear();
	return position;
}
}