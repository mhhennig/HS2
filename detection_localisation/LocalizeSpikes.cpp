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
	centered_amps: deque<tuple<int, int>>
		A deque containing all non-zero amplitudes and their neighbors. Used for center of mass.

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
		channel = get<0>(centered_amps.at(i));
		X_coordinate = Parameters::channel_positions[channel][0];
		Y_coordinate = Parameters::channel_positions[channel][1];
		if(curr_amp < 0) {
            weight = 0;
        } else {
            weight = curr_amp;
        }
		X_numerator += weight * X_coordinate;
		Y_numerator += weight * Y_coordinate;
		denominator += weight;
	}

	X = (float)(X_numerator) / (float)(denominator);
	Y = (float)(Y_numerator) / (float)(denominator);

    if(X <= 0 || Y <= 0) {
        for(int i = 0; i < centered_amps_size; i++) {
    		curr_amp = get<1>(centered_amps.at(i));
    		channel = get<0>(centered_amps.at(i));
    		X_coordinate = Parameters::channel_positions[channel][0];
    		Y_coordinate = Parameters::channel_positions[channel][1];
            weight = curr_amp;
            cout << "X coord: " << X_coordinate << endl;
            cout << "Y coord: " << Y_coordinate << endl;
            cout << "Weight: " << weight << endl;
    	}
        cout << "X: " << X << endl;
        cout << "Y: " << Y << endl;
     }

	return make_tuple(X, Y);
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
	deque<tuple<int, int>> amps;
	int curr_largest_amp = INT_MIN; //arbitrarily small to make sure that it is immediately overwritten
	int curr_neighbor_channel;
	int curr_amp;

	int amp_cutout_size = spike_to_be_localized.amp_cutouts.size();
	int neighbor_frame_span = Parameters::spike_delay * 2 + 1;
	int neighbor_count = amp_cutout_size / neighbor_frame_span;
	int matrix_offset = 0;

	for (int i = 0; i < neighbor_count; i++) {
		curr_neighbor_channel = Parameters::inner_neighbor_matrix[spike_to_be_localized.channel][i];
        if(Parameters::masked_channels[curr_neighbor_channel] != 0) {
    		for (int j = 0; j < neighbor_frame_span; j++) {
    			curr_amp = spike_to_be_localized.amp_cutouts.at(j + matrix_offset);
                if(curr_amp > 100000) {
                    cout << "CURR AMP TOO BIG: " << curr_amp << endl;
                }
    			if(curr_amp > curr_largest_amp) {
    				curr_largest_amp = curr_amp;
    			}
    		}
            // if(curr_largest_amp < 0) {
            //     for (int j = 0; j < neighbor_frame_span; j++) {
        	// 		curr_amp = spike_to_be_localized.amp_cutouts.at(j + matrix_offset);
        	// 		cout << curr_amp << endl;
        	// 	}
            // }
    		amps.push_back(make_tuple(curr_neighbor_channel, curr_largest_amp));
    		curr_largest_amp = INT_MIN;
    		matrix_offset += neighbor_frame_span;
        }
	}
    // int do_median = 0;
    // int median = 0;
    // int amps_size = amps.size();
    // if(do_median == 1) {
    //     sort(begin(amps), end(amps), CustomLessThan()); //sort the array
    //     //Find median of array
    // 	if(amps_size % 2 == 0) {
    // 		median = (get<1>(amps.at(amps_size/2)) + get<1>(amps.at(amps_size/2 + 1)))/2;
    // 	}
    // 	else {
    // 		median = get<1>(amps.at(amps_size/2));
    // 	}
    // }
	// //Center amplitudes
    // deque<tuple<int, int>> centered_amps;
    // if(amps_size != 1) {
    // 	for(int i = 0; i < amps_size; i++) {
    // 		int curr_neighbor = get<0>(amps.at(i));
    // 		int curr_amp = get<1>(amps.at(i));
    // 		int new_amp = curr_amp - median;
    // 		//if(new_amp > 0) {
    // 			//centered_amps.push_back(make_tuple(curr_neighbor, new_amp));
    // 		//}
    //         centered_amps.push_back(make_tuple(curr_neighbor, new_amp));
    // 	}
    // } else {
    //     centered_amps.push_back(amps.at(0));
    // }
    //
    // int centered_amps_size = centered_amps.size();
    // //Subtracting median made all values equal to or smaller than 0
    // if(centered_amps_size == 0) {
    //     centered_amps = amps;
    // }
	tuple<float,float> position = centerOfMass(amps);
	amps.clear();
	//centered_amps.clear();
	return position;
}

}
