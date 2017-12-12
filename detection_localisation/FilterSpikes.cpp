#include "FilterSpikes.h"

using namespace std;

namespace FilterSpikes {

Spike filterSpikes(Spike first_spike, ofstream& filteredsp) {
	/*Removes all duplicate spikes and returns the max spike

	Parameters
	----------
	first_spike: Spike
		The first spike detected in the event.

    Returns
    -------
    max_spike: Spike
        The largest amplitude spike belonging to the event of the first spike.

    */
    Spike max_spike;
    max_spike = findMaxSpikeNeighbor(first_spike);
    if(Parameters::spikes_to_be_processed.size() != 0) {
        max_spike = updateNeighborsMaxSpike(max_spike);
        //Find the max amplitude neighbor of the first spike
        filterOuterNeighbors(max_spike, filteredsp);
        filterInnerNeighbors(max_spike, filteredsp);
    }
    return max_spike;
}

Spike findMaxSpikeNeighbor(Spike first_spike) {
    /*Finds the max amplitude neighbor of the first spike in the current window

	Parameters
	----------
	first_spike: Spike
		The first spike detected in the event.

    Returns
    -------
	max_spike: Spike
		The largest amplitude spike belonging to the event of the first spike.
	*/

    int curr_channel, curr_amp, curr_frame;
    Spike curr_spike;
    Spike max_spike;
    int max_spike_amp;
    max_spike = first_spike;
    max_spike_amp = first_spike.amplitude;
    //Find the max amplitude neighbor of the first spike
    deque<Spike>::iterator it;
    deque<Spike>::iterator it_final;
	it = Parameters::spikes_to_be_processed.begin();
    int index = it - Parameters::spikes_to_be_processed.begin();

    //Find max
	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		curr_amp = it->amplitude;
		curr_frame = it->frame;
		if(areNeighbors(first_spike.channel, curr_channel)) {
            if(curr_amp >= max_spike_amp) {
                if(curr_frame <=  first_spike.frame + Parameters::noise_duration) {
                    max_spike = curr_spike;
                    max_spike_amp = curr_amp;
                    index = it - Parameters::spikes_to_be_processed.begin();
                }
            }
	    }
        ++it;
	}
    deque<Spike>::iterator max_it = Parameters::spikes_to_be_processed.begin() + index;
    Parameters::spikes_to_be_processed.erase(max_it);
    return max_spike;
}

Spike updateNeighborsMaxSpike(Spike max_spike) {
    /*Finds and updates all the inner and outer neighbors that spiked of the max amplitude
    spike.

	Parameters
	----------
	max_spike: Spike
        The largest amplitude spike belonging to the event of the first spike.

    Returns
    -------
    max_spike: Spike
        The largest amplitude spike belonging to the event of the first spike.
        Now its inner and outer neighbors are updated with spike information.
	*/

    int curr_channel, curr_amp, curr_frame;
    Spike curr_spike;
    vector<int> spiking_neighbors;

    //Fill all neighbors that spiked with amplitude of their spike
    deque<Spike>::iterator it;
	it = Parameters::spikes_to_be_processed.begin();
	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		curr_amp = it->amplitude;
        curr_frame = it->frame;
        if(areNeighbors(max_spike.channel, curr_channel)) {
            if(isInnerNeighbor(max_spike.channel, curr_channel)) {
                Event event;
        		event.channel = curr_channel;
        		event.frame = curr_frame;
        		event.amplitude = curr_amp;
                max_spike.inner_neighbors.push_back(event);
                spiking_neighbors.push_back(curr_channel);
            }
            else {
                Event event;
        		event.channel = curr_channel;
        		event.frame = curr_frame;
        		event.amplitude = curr_amp;
                max_spike.outer_neighbors.push_back(event);
                spiking_neighbors.push_back(curr_channel);
            }
	    }
        ++it;
	}

    //Fill all neighbors that didn't spike with zero amplitude
    int curr_neighbor;
    for(int i = 0; i < Parameters::max_neighbors; i++) {
        curr_neighbor = Parameters::neighbor_matrix[max_spike.channel][i];
        if(std::find(spiking_neighbors.begin(), spiking_neighbors.end(), curr_neighbor) != spiking_neighbors.end()) {
            //Already updated
        }
        else {
            if(max_spike.channel != curr_neighbor) {
                if(isInnerNeighbor(max_spike.channel, curr_neighbor)) {
                    Event event;
            		event.channel = curr_neighbor;
            		event.frame = -1;
            		event.amplitude = 0;
                    max_spike.inner_neighbors.push_back(event);
                }
                else {
                    Event event;
            		event.channel = curr_neighbor;
            		event.frame = -1;
            		event.amplitude = 0;
                    max_spike.outer_neighbors.push_back(event);
                }
            }
        }
    }

    return max_spike;
}

void filterOuterNeighbors(Spike max_spike, ofstream& filteredsp) {
    /*Filters or leaves all outer neighbors based on minimum spanning tree
    algorithm. Basically, the outer spikes must pass through an inner spike
    to reach the maximum

	Parameters
	----------
	max_spike: Spike
		The max spike detected in the event.
	*/

    Spike curr_spike;
	int curr_channel, curr_frame;
    vector<Spike> outer_spikes_to_be_filtered;
    deque<Spike>::iterator it;
	it = Parameters::spikes_to_be_processed.begin();

	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = curr_spike.channel;
		if(areNeighbors(max_spike.channel, curr_channel)) {
            if(!isInnerNeighbor(max_spike.channel, curr_channel)) {
                if(filteredOuterSpike(curr_spike, max_spike)) {
                    //filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " PN ratio: " << posToNegRatio(curr_spike) << " Area/Amp: " << areaUnderSpike(curr_spike) << " RP time: " << repolarizationTime(curr_spike) << " Filtered by " << max_spike.channel << endl;
                    filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << "  " << posToNegRatio(curr_spike)  << " " << areaUnderSpike(curr_spike) << " " << repolarizationTime(curr_spike) << endl;
                    outer_spikes_to_be_filtered.push_back(curr_spike);
                    ++it;
                }
                else {
                    //Not a decaying outer spike (probably a new spike), filter later
                    ++it;
                }
            }
            else {
                //Inner neighbor - all inner neighbors that need filtering will be filtered later
                ++it;
            }
		}
		else {
            //Not a neighbor, filter later
			++it;
		}
	}

    //Filter all spikes that need to be filtered all together at once
    Spike curr_spike_to_be_filtered;
    int curr_channel_to_be_filtered;
    int curr_frame_to_be_filtered;
    vector<Spike>::iterator it2;
    it2 = outer_spikes_to_be_filtered.begin();
    while(it2 != outer_spikes_to_be_filtered.end())
	{
        curr_spike_to_be_filtered = *it2;
        curr_channel_to_be_filtered = curr_spike_to_be_filtered.channel;
        curr_frame_to_be_filtered = curr_spike_to_be_filtered.frame;

    	it = Parameters::spikes_to_be_processed.begin();
        while(it != Parameters::spikes_to_be_processed.end())
    	{
            curr_spike = *it;
    		curr_channel = curr_spike.channel;
            curr_frame = curr_spike.frame;
            if(curr_channel == curr_channel_to_be_filtered && curr_frame == curr_frame_to_be_filtered) {
                it = Parameters::spikes_to_be_processed.erase(it);
                break;
            }
            else {
                ++it;
            }
        }
        ++it2;
	}
}

bool filteredOuterSpike(Spike outer_spike, Spike max_spike) {
    //True if should be filtered, false by default
    bool filtered_spike = false;
    bool IS_INNER_EVENT = true;
    bool shares_inner = false;
    int NOT_A_CHANNEL = -10;

    //Checks to see if outer_spike shares an inner neighbor with max_spike (Base Case)
    int curr_inner_neighbor;
    for(int i = 0; i < Parameters::max_neighbors - 1; i++) {
        curr_inner_neighbor = Parameters::inner_neighbor_matrix[outer_spike.channel][i];
        if(curr_inner_neighbor == -1) {
            break;
        }
        else {
            Event curr_max_inner_event = getEventfromChannel(curr_inner_neighbor, max_spike, IS_INNER_EVENT);
            if(curr_max_inner_event.channel != NOT_A_CHANNEL) {
                if(outer_spike.amplitude < curr_max_inner_event.amplitude*Parameters::noise_amp_percent) {
                    if(outer_spike.frame < curr_max_inner_event.frame - Parameters::noise_duration) {
                        //outer spike occurs too far before inner spike, probably new spike
                        shares_inner = true;
                    }
                    else {
                        //Shares an inner neighbor, spikes at a reasonable time, amplitude has decayed enough over this distance, filter
                        filtered_spike = true;
                        shares_inner = true;
                        break;
                    }
                }
                else {
                    //amplitude too big to be a duplicate spike
                    shares_inner = true;
                }
            }
            else {
                //Keep looking for shared inner neighbor
            }
        }
    }
    if(shares_inner == true) {
        return filtered_spike;
    }
    else {
        //Doesn't share an inner_neighbor with max_spike (find closest inner neighbor to centeral channel and its amplitude)
        int closest_inner_neighbor_channel = getClosestInnerNeighborChannel(outer_spike.channel, max_spike.channel);
        Event closest_event = getEventfromChannel(closest_inner_neighbor_channel, max_spike, !IS_INNER_EVENT);
        int closest_amp = closest_event.amplitude;
        //Outer spike amplitude too big, probably new spike, filter later
        if(outer_spike.amplitude >= closest_amp*Parameters::noise_amp_percent) {
            return filtered_spike;
        }
        else {
            //Find the closest neighbor in the spike vector and recurse with it (Recursive Step)
            Spike closest_spike = getSpikefromChannel(closest_inner_neighbor_channel);
            if(outer_spike.frame >=  closest_spike.frame - Parameters::noise_duration) {
                return filteredOuterSpike(closest_spike, max_spike);
            }
            else {
                //Occured too far before the inner spike to be related
                return filtered_spike;
            }
        }
        //Should never come to this
        return filtered_spike;
    }
}

int getClosestInnerNeighborChannel(int outer_channel, int central_channel) {
    int closest_dist = 10000;
    int closest_inner_neighbor_channel = -1;
    int curr_dist;
    int curr_inner_channel;
    for(int i = 0; i < Parameters::max_neighbors - 1; i++) {
        curr_inner_channel = Parameters::inner_neighbor_matrix[outer_channel][i];
        if(curr_inner_channel == -1) {
            break;
        } else {
            curr_dist = channelsDist(curr_inner_channel, central_channel);
            if(curr_dist < closest_dist) {
                closest_dist = curr_dist;
                closest_inner_neighbor_channel = curr_inner_channel;
            }
        }
    }
    return closest_inner_neighbor_channel;
}

Spike getSpikefromChannel(int spike_channel) {
    int curr_channel;
    Spike curr_spike;
    Spike channel_spike;
    deque<Spike>::iterator it;
    it = Parameters::spikes_to_be_processed.begin();
    while(it != Parameters::spikes_to_be_processed.end())
    {
        curr_spike = *it;
        curr_channel = it->channel;
        if(curr_channel == spike_channel) {
            channel_spike = curr_spike;
            break;
        }
        else {
            ++it;
        }
    }
    return channel_spike;
}

Event getEventfromChannel(int channel, Spike curr_spike, bool is_inner_event) {
    Event event_channel;
    //If returned channel has -10 as channel number, then no event was found.
    event_channel.channel = -10;
    vector<Event>::iterator it;
    vector<Event>::iterator end;
    Event curr_event;
    int curr_channel_event;
    if(is_inner_event) {
        it = curr_spike.inner_neighbors.begin();
        end = curr_spike.inner_neighbors.end();
    }
    else {
        it = curr_spike.outer_neighbors.begin();
        end = curr_spike.outer_neighbors.end();
    }
    while(it != end)
    {
        curr_event = *it;
        curr_channel_event = curr_event.channel;
        if(channel == curr_channel_event) {
            event_channel = curr_event;
            break;
        }
        else {
            ++it;
        }
    }
    return event_channel;
}

void filterInnerNeighbors(Spike max_spike, ofstream& filteredsp) {
    /*Filters or leaves all outer neighbors based on minimum spanning tree
    algorithm. Basically, the outer spikes must pass through an inner spike
    to reach the maximum

	Parameters
	----------
	max_spike: Spike
		The max spike detected in the event.
	*/

	Spike curr_spike;
	int curr_channel, curr_amp, curr_frame;
    deque<Spike>::iterator it;
	it = Parameters::spikes_to_be_processed.begin();

	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		curr_amp = it->amplitude;
        curr_frame = it->frame;
		if(areNeighbors(max_spike.channel, curr_channel)) {
            if(isInnerNeighbor(max_spike.channel, curr_channel)) {
                if(curr_amp < max_spike.amplitude) {
                    if(curr_amp >=  max_spike.amplitude*Parameters::noise_amp_percent) {
                        if(curr_frame >  max_spike.frame + Parameters::noise_duration) {
                            //Inner spike occurred much later that max_spike and is only slightly smaller (probably new spike), filter later
                            ++it;
                        }
                        else {
                            //filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " PN ratio: " << posToNegRatio(curr_spike)  << " Area/Amp:: " << areaUnderSpike(curr_spike) << " RP time: " << repolarizationTime(curr_spike) << " Filtered by " << max_spike.channel << endl;
                            filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << "  " << posToNegRatio(curr_spike)  << " " << areaUnderSpike(curr_spike) << " " << repolarizationTime(curr_spike) << endl;
                            it = Parameters::spikes_to_be_processed.erase(it);
                        }
                    }
                    else {
                        //filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " PN ratio: " << posToNegRatio(curr_spike) << " Area/Amp:: " << areaUnderSpike(curr_spike) << " RP time: " << repolarizationTime(curr_spike) <<  " Filtered by " << max_spike.channel << endl;
                        filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << "  " << posToNegRatio(curr_spike)  << " " << areaUnderSpike(curr_spike) << " " << repolarizationTime(curr_spike) << endl;
                        it = Parameters::spikes_to_be_processed.erase(it);
                    }
                }
                else {
                    //Inner neighbor has larger amplitude that max neighbor (probably new spike), filter later
                    ++it;
                }
            }
            else {
                //In outer neighbors, filter later
                ++it;
            }
		}
		else {
            //Not a neighbor, filter later
			++it;
		}
	}
}

float channelsDist(int start_channel, int end_channel) {
    /*Finds the distance between two channels

    Parameters
    ----------
    start_channel: int
        The start channel where distance measurement begins
    end_channel: int
        The end channel where distance measurement ends

    Returns
    -------
    dist: float
        The distance between the two channels
    */
    float start_position_x;
    float start_position_y;
    float end_position_x;
    float end_position_y;
    float x_displacement;
    float y_displacement;
    float dist;

    start_position_x = Parameters::channel_positions[start_channel][0];
    start_position_y = Parameters::channel_positions[start_channel][1];
    end_position_x = Parameters::channel_positions[end_channel][0];
    end_position_y = Parameters::channel_positions[end_channel][1];
    x_displacement = start_position_x - end_position_x;
    y_displacement = start_position_y - end_position_y;
    dist = sqrt(pow(x_displacement, 2) + pow(y_displacement, 2));

    return dist;
}

bool isInnerNeighbor(int center_channel, int curr_channel) {
    /*Determines whether or not two channels are inner neighbors.

    Parameters
    ----------
    center_channel: int
        The channel number whose inner neighborhoold is checked.
    curr_channel: int
        The channel number of the potential neighbor.

    Returns
    -------
    is_inner_neighbor: bool
        True if the channels are inner neighbors.
        False if the channels are not inner neighbors.
    */
    bool is_inner_neighbor = false;

    int curr_inner_neighbor;
    for(int i = 0; i < Parameters::max_neighbors; i++) {
        curr_inner_neighbor = Parameters::inner_neighbor_matrix[center_channel][i];
        if(curr_inner_neighbor == curr_channel) {
            is_inner_neighbor = true;
            break;
        }
        else if(curr_inner_neighbor == -1) {
            break;
        }
    }
    return is_inner_neighbor;
}

bool areNeighbors(int channel_one, int channel_two) {
	/*Determines whether or not two channels are neighbors.

	Parameters
	----------
	channel_one: int
		The channel number whose neighborhoold is checked.
	channel_two: int
		The channel number of the potential neighbor.

	Returns
	-------
	are_neighbors: bool
		True if the channels are neighbors.
		False if the channels are not neighbors.
	*/
	bool are_neighbors = false;
    int curr_neighbor;
	for(int i = 0; i < Parameters::max_neighbors; i++) {
        curr_neighbor = Parameters::neighbor_matrix[channel_one][i];
		if(curr_neighbor == channel_two) {
			are_neighbors = true;
			break;
		}
        else if(curr_neighbor == -1) {
            break;
        }
	}
	return are_neighbors;
}

float posToNegRatio(Spike spike) {
   int32_t most_neg_reading = 2147483647;
   int32_t most_pos_reading = -2147483647;
   int32_t curr_written_reading;
   int index_neg_spike = 0;
   int index_spike = 0;
   for(int i = 0; i < Parameters::cutout_start + Parameters::noise_duration + 1; i++) {
       curr_written_reading = spike.written_cutout.at(i);;
       if(curr_written_reading < most_neg_reading) {
           most_neg_reading = curr_written_reading;
           index_neg_spike = index_spike;
       }
       ++index_spike;
    }


   for(int i = 0; i < index_neg_spike; i++) {
       curr_written_reading = spike.written_cutout.at(i);
       if(curr_written_reading > most_pos_reading) {
           most_pos_reading = curr_written_reading;
       }
    }

   float pos_to_neg_ratio = -float(most_pos_reading)/float(most_neg_reading);
   return pos_to_neg_ratio;
}

float repolarizationTime(Spike spike) {
   int32_t most_neg_reading = 2147483647;
   int32_t most_pos_reading = -2147483647;
   int32_t curr_written_reading;
   int32_t prev_written_reading;
   int index_neg_spike = 0;
   int index_spike = 0;
   for(int i = 0; i < Parameters::cutout_start + Parameters::noise_duration + 1; i++) {
       curr_written_reading = spike.written_cutout.at(i);;
       if(curr_written_reading < most_neg_reading) {
           most_neg_reading = curr_written_reading;
           index_neg_spike = index_spike;
       }
       ++index_spike;
    }

   int frames = 0;
   for(size_t i = index_neg_spike; i < spike.written_cutout.size(); i++) {
       curr_written_reading = spike.written_cutout.at(i);
       if(i == index_neg_spike) {
           prev_written_reading = curr_written_reading;
           ++frames;
       }
       else {
           if(curr_written_reading > 0 && curr_written_reading < prev_written_reading) {
               break;
           }
           else {
               prev_written_reading = curr_written_reading;
               ++frames;
           }
       }
    }

   return frames - 1;
}

double areaUnderSpike(Spike spike) {
    int32_t most_neg_reading = 2147483647;
    int32_t curr_written_reading;
    int index_neg_spike = 0;
    int index_spike = 0;
    for(int i = 0; i < Parameters::cutout_start + Parameters::noise_duration + 1; i++) {
        curr_written_reading = spike.written_cutout.at(i);;
        if(curr_written_reading < most_neg_reading) {
            most_neg_reading = curr_written_reading;
            index_neg_spike = index_spike;
        }
        ++index_spike;
     }

    double sum = 0.0, trapz;
    int h = 1;
    int amp_cutout_size = Parameters::cutout_start*2 + 1;
    for (int i = 0; i < amp_cutout_size; ++i){
        if ( i == 0 || i == amp_cutout_size-1 ) {
            sum += ((double) spike.written_cutout.at(i))/2.0;
        }// for the first and last elements
        else {
            sum += ((double) spike.written_cutout.at(i)); // the rest of data
        }
    }
    trapz = sum*h; // the result

    if(most_neg_reading == 0) {
        return trapz;
    }
    else {
        return trapz/most_neg_reading;
    }
}


}
