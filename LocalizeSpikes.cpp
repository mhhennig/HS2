#include "LocalizeSpikes.h"

tuple<float,float> centerOfMass(vector<float>* amps, int spike_channel)
{
	float X = 0;
	float Y = 0;
	float X_numerator = 0;
	float Y_numerator = 0;
	float denominator = 0;
	int X_coordinate;
	int Y_coordinate;
	int channel;
	float weight;
	int amps_size = amps->size();

	for(int i = 0; i < amps_size; i++) {
		if(amps->at(i) < 0) {
			continue;
		}
		else {
			channel = Parameters::neighbor_matrix[spike_channel][i];
			X_coordinate = Parameters::channel_positions[channel][0];
			Y_coordinate = Parameters::channel_positions[channel][1];
			weight = amps->at(i);
			X_numerator += weight * X_coordinate;
			Y_numerator += weight * Y_coordinate;
			denominator += weight;
		}
	}

	X = X_numerator / denominator;
	Y = Y_numerator / denominator;

	return make_tuple(X, Y);
}

int getNumNeighbors(int channel)
{
	int num_neighbors = 0;
	for(int i = 0; i < Parameters::max_neighbors; i++) {
		if(Parameters::neighbor_matrix[channel][i] != -1) {
			num_neighbors += 1;
		}
		else {
			break;
		}
	}
	return num_neighbors;
}

tuple<float, float> localizeSpike(Spike spike_to_be_localized)
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
	int num_neighbors = getNumNeighbors(spike_channel);
	vector<float> amps;
	float curr_largest_amp = -100000; //arbitrarily small to make sure that it is immediately overwritten
	int curr_reading, start_cutout, end_cutout, curr_neighbor_channel;
	float curr_amp;

	for(int i = 0; i < num_neighbors; i++) {
		curr_neighbor_channel = Parameters::neighbor_matrix[spike_channel][i];
		start_cutout = spike_start_frame*Parameters::num_channels + curr_neighbor_channel;
		end_cutout = spike_end_frame*Parameters::num_channels + curr_neighbor_channel;
		for(int j = 0; j < end_cutout; j++) {
			curr_reading = Parameters::raw_data[start_cutout + j];
			curr_amp = float(-1*(curr_reading - Parameters::aGlobal * ASCALE - Parameters::baselines[curr_neighbor_channel]));
			if(curr_amp > curr_largest_amp) {
				curr_largest_amp = curr_amp;
			}
		}
		amps.push_back(curr_largest_amp);
	}

	sort(begin(amps), end(amps)); //sort the array
	
	//Find median of array
	int amps_size = amps.size();
	float median;
	if(amps_size % 2 == 0) {
		median = (amps.at(amps_size/2) + amps.at(amps_size/2 + 1))/2.0;
	}
	else {
		median = amps.at(amps_size/2);	
	}

	//Center amplitudes
	for(int i = 0; i < amps_size; i++) {
		amps.at(i) = amps.at(i) - median;
	}

	tuple<float,float> position = centerOfMass(&amps, spike_channel);

	return position;
}