#include "LocalizeSpikes.h"

namespace LocalizeSpikes {

struct CustomLessThan {
  bool operator()(tuple<int, int> const &lhs,
                  tuple<int, int> const &rhs) const {
    return std::get<1>(lhs) < std::get<1>(rhs);
  }
};

tuple<float, float> localizeSpike(Spike spike_to_be_localized) {
  /*Estimates the X and Y position of where a spike occured on the probe.

     Parameters
     ----------
     spike_to_be_localized: Spike
     The spike that will be used to determine where the origin of the spike
     occurred.

     Returns
     -------
     position: tuple<float, float>
     An X and Y coordinate tuple that corresponds to where the spike occurred.
   */

  vector<int> waveforms = get<0>(spike_to_be_localized.waveformscounts);
  deque<tuple<tuple<float, float>, int>> com_positions_amps;
  int matrix_offset = 0;
  int curr_largest_amp;
  int curr_neighbor_channel;
  int curr_amp;
  int sum_amp;
  int neighbor_count;
  int cutout_size;
  int curr_max_channel;
  for (int i = 0; i < Parameters::num_com_centers; i++) {
    deque<tuple<int, int>> amps;
    curr_largest_amp = INT_MIN; // arbitrarily small to make sure that it is
                                    // immediately overwritten
    neighbor_count = get<1>(spike_to_be_localized.waveformscounts)[i];
    cutout_size = Parameters::noise_duration*2;
    curr_max_channel = spike_to_be_localized.largest_channels[i];
    // compute amplitudes using sum over 2*noise_duration data points
    for (int j = 0; j < neighbor_count; j++) {
      curr_neighbor_channel =
          Parameters::inner_neighbor_matrix[curr_max_channel][j];
      if (Parameters::masked_channels[curr_neighbor_channel] == 1) {
	sum_amp = 0;
       for (int k = 0; k < cutout_size; k++) {
          sum_amp += waveforms[k + matrix_offset];
/*
          curr_amp = waveforms[k + matrix_offset];
          if (curr_amp > curr_largest_amp) {
            sum_amp += curr_amp;
            curr_largest_amp = curr_amp;
          }
*/
        }
	//curr_largest_amp = curr_amp;
        amps.push_back(make_tuple(curr_neighbor_channel, sum_amp));
        //curr_largest_amp = INT_MIN;
        matrix_offset += cutout_size;
      }
    }

    // compute median, threshold at median
    int do_correction = 1;
    int correct = 0;
    int amps_size = amps.size();
    if (Parameters::debug) {
      cout << "Correction phase..." << endl;
    }
    if (do_correction == 1) {
      sort(begin(amps), end(amps), CustomLessThan()); // sort the array
      //correct = get<1>(amps.at(0))-1;
      if (Parameters::debug) {
              cout << "Amps size: " << amps_size << endl;
      }
      if(amps_size % 2 == 0) {
              correct = (get<1>(amps.at(amps_size/2 - 1)) + get<1>(amps.at(amps_size/2)))/2;
      }
      else {
              correct = get<1>(amps.at(amps_size/2));
      }

    }
    if (Parameters::debug) {
      cout << "Done correcting..." << endl;
    }
    // Correct amplitudes (threshold)
    deque<tuple<int, int>> centered_amps;
    if (amps_size != 1) {
      for (int i = 0; i < amps_size; i++) {
	if(get<1>(amps.at(i)) - correct>0) {
        	centered_amps.push_back(
            	make_tuple(get<0>(amps.at(i)), get<1>(amps.at(i)) - correct));
      		}
	}
    } else {
      centered_amps.push_back(amps.at(0));
    }

    tuple<float, float> position = centerOfMass(centered_amps);
    tuple<tuple<float, float>, int> position_amp_tuple = make_tuple(position, 1);
    com_positions_amps.push_back(position_amp_tuple);

    amps.clear();
    centered_amps.clear();
  }

  tuple<float, float> reweighted_com;
  if(com_positions_amps.size() > 1) {
    reweighted_com = reweightedCenterOfMass(com_positions_amps);
  } else {
    reweighted_com = get<0>(com_positions_amps[0]);
  }

  return reweighted_com;
}

tuple<float, float> reweightedCenterOfMass(deque<tuple<tuple<float, float>, int>> com_positions_amps) {
  float X = 0;
  float Y = 0;
  float X_numerator = 0;
  float Y_numerator = 0;
  int denominator = 0;
  float X_coordinate;
  float Y_coordinate;
  int weight; // contains the amplitudes for the center of mass calculation.
              // Updated each localization

  for (int i = 0; i < Parameters::num_com_centers; i++) {
    X_coordinate = get<0>(get<0>((com_positions_amps[i])));
    Y_coordinate = get<1>(get<0>((com_positions_amps[i])));
    weight = get<1>(com_positions_amps[i]);
    if (weight < 0) {
      cout << "\ncenterOfMass::weight < 0 - this should not happen\n";
    }
    X_numerator += weight * X_coordinate;
    Y_numerator += weight * Y_coordinate;
    denominator += weight;
  }

  if(denominator == 0) {
    cout << "Whopodis" << endl;
    for (int i = 0; i < Parameters::num_com_centers; i++) {
      X_coordinate = get<0>(get<0>((com_positions_amps[i])));
      Y_coordinate = get<1>(get<0>((com_positions_amps[i])));
      weight = get<1>(com_positions_amps[i]);
      if (weight < 0) {
        cout << "\ncenterOfMass::weight < 0 - this should not happen\n";
      }
      cout << "Weight" << weight << endl;
      cout << "X coordinate" << X_coordinate << endl;
      cout << "Y coordinate" << Y_coordinate << endl;
    }
  }

  X = (X_numerator) / (float)(denominator);
  Y = (Y_numerator) / (float)(denominator);

  return make_tuple(X, Y);
}

tuple<float, float> centerOfMass(deque<tuple<int, int>> centered_amps) {
  /*Calculates the center of mass of a spike to calculate where it occurred
     using a weighted average.

     Parameters
     ----------
     centered_amps: deque<tuple<int, int>>
     A deque containing all non-zero amplitudes and their neighbors. Used for
     center of mass.

     Returns
     -------
     position: tuple<float, float>
     An X and Y coordinate tuple that corresponds to where the spike occurred.
   */
  // int curr_amp;
  float X = 0;
  float Y = 0;
  float X_numerator = 0;
  float Y_numerator = 0;
  int denominator = 0;
  float X_coordinate;
  float Y_coordinate;
  int channel;
  int weight; // contains the amplitudes for the center of mass calculation.
              // Updated each localization
  int centered_amps_size = centered_amps.size();

  if (Parameters::debug) {
    cout << "Done localizing..." << endl;
  }
  for (int i = 0; i < centered_amps_size; i++) {
    weight = get<1>(centered_amps.at(i));
    channel = get<0>(centered_amps.at(i));
    X_coordinate = Parameters::channel_positions[channel][0];
    Y_coordinate = Parameters::channel_positions[channel][1];
    if (weight < 0) {
      cout << "\ncenterOfMass::weight < 0 - this should not happen\n";
    }
    X_numerator += weight * X_coordinate;
    Y_numerator += weight * Y_coordinate;
    denominator += weight;
  }

  if (denominator == 0) //| (X>10 & X<11))
  {
    // cout << "\ncenterOfMass::denominator == 0 - This should not happen\n";
    for (int i = 0; i < centered_amps_size; i++) {
      channel = get<0>(centered_amps.at(i));
      // cout << " " << get<1>(centered_amps.at(i)) << " "
      //      << Parameters::channel_positions[channel][0] << "\n";
      X = Parameters::channel_positions[channel][0];
      Y = Parameters::channel_positions[channel][1];
    }
    cout << "\n";
  }
  else {
    X = (float)(X_numerator) / (float)(denominator);
    Y = (float)(Y_numerator) / (float)(denominator);
  }

  return make_tuple(X, Y);
}

}
