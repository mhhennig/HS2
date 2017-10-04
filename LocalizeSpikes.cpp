#include "LocalizeSpikes.h"

namespace LocalizeSpikes {

struct CustomLessThan
{
    bool operator()(tuple<int, int> const &lhs, tuple<int, int> const &rhs) const
    {
        return std::get<1>(lhs) < std::get<1>(rhs);
    }
};

tuple<float, float> centerOfMass(deque<tuple<int, int>> centered_amps)
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
	int curr_amp;
	float X = 0;
	float Y = 0;
	int X_numerator = 0;
	int Y_numerator = 0;
	int denominator = 0;
	int X_coordinate;
	int Y_coordinate;
	int channel;
	int weight; //contains the amplitudes for the center of mass calculation. Updated each localization
	int centered_amps_size = centered_amps.size();

	for(int i = 0; i < centered_amps_size; i++) {
		curr_amp = get<1>(centered_amps.at(i));
		if(curr_amp > 0) {
			channel = get<0>(centered_amps.at(i));
			X_coordinate = Parameters::channel_positions[channel][0];
			Y_coordinate = Parameters::channel_positions[channel][1];
			weight = curr_amp;
			//cout << "Channel: " << channel << endl;
			//cout << "X_coordinate: " << X_coordinate << endl;
			//cout << "Y_coordinate: " << Y_coordinate<< endl;
			//cout << "weight:  " << curr_amp << endl;
			X_numerator += weight * X_coordinate;
			Y_numerator += weight * Y_coordinate;
			denominator += weight;
		}
	}

	X = (float)(X_numerator) / (float)(denominator);
	Y = (float)(Y_numerator) / (float)(denominator);

	return make_tuple(X, Y);
}

tuple<float, float> localizeSpike(Spike spike_to_be_localized)
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
	deque<tuple<int, int>> amps;
	int spike_channel = spike_to_be_localized.channel;
	int curr_largest_amp = -100000; //arbitrarily small to make sure that it is immediately overwritten
	int curr_reading, curr_neighbor_channel;
	int curr_amp;
	//cout << "Localizing Spike: " << spike_to_be_localized.channel << " " << spike_to_be_localized.frame << " " << spike_to_be_localized.amplitude << '\n';

	int amp_cutout_size = spike_to_be_localized.amp_cutouts.size();
	for(int i = 0; i < amp_cutout_size; i++) {
		curr_neighbor_channel = Parameters::neighbor_matrix[spike_channel][i /(Parameters::spike_delay*2 + 1)];
		curr_amp = spike_to_be_localized.amp_cutouts.at(i);
		//cout << curr_amp << endl;
		//cout << "curr reading: " << curr_reading << endl;
		//cout << "a global: " << Parameters::aGlobal << endl;
		//cout << "baseline: " << Parameters::baselines[curr_neighbor_channel][(Parameters::index_baselines - baseline_frame)] << endl;
		//curr_amp = -1*(curr_reading - Parameters::aGlobal * ASCALE - Parameters::baselines[curr_neighbor_channel][(Parameters::index_baselines - baseline_frame) % (Parameters::spike_delay + 1)]);
		if(curr_amp > curr_largest_amp) {
			curr_largest_amp = curr_amp;
		}
		if(i % (Parameters::spike_delay*2 + 1) == Parameters::spike_delay*2) {
			//cout << "curr neighbor: " << curr_neighbor_channel << " largest amp: " << curr_largest_amp << endl;
			amps.push_back(make_tuple(curr_neighbor_channel, curr_largest_amp));
			curr_largest_amp = -10000;
		}
	}

	sort(begin(amps), end(amps), CustomLessThan()); //sort the array
	
	//Find median of array
	int amps_size = amps.size();
	int median;
	if(amps_size % 2 == 0) {
		median = (get<1>(amps.at(amps_size/2)) + get<1>(amps.at(amps_size/2 + 1)))/2;
	}
	else {
		median = get<1>(amps.at(amps_size/2));	
	}

	//Center amplitudes
	deque<tuple<int, int>> centered_amps;
	for(int i = 0; i < amps_size; i++) {
		int curr_neighbor = get<0>(amps.at(i));
		int curr_amp = get<1>(amps.at(i));
		int new_amp = curr_amp - median;
		centered_amps.push_back(make_tuple(curr_neighbor, new_amp));
		//cout << "centered amp: " << get<1>(centered_amps.at(i)) << endl;
	}
	
	tuple<float,float> position = centerOfMass(centered_amps);
	amps.clear();
	centered_amps.clear();
	return position;
}

}