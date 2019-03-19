#include "FilterSpikes.h"

using namespace std;

namespace FilterSpikes {

Spike filterSpikesDecay(Spike first_spike, ofstream& filteredsp) {
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
    //Find the max amplitude neighbor of the first spike
    max_spike = findMaxSpikeNeighbor(first_spike);
    if(Parameters::spikes_to_be_processed.size() != 0) {
        filterOuterNeighbors(max_spike, filteredsp);
        filterInnerNeighbors(max_spike, filteredsp);
    }
    return max_spike;
}

Spike filterSpikesAll(Spike first_spike, ofstream& filteredsp) {
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
    //Find the max amplitude neighbor of the first spike
    max_spike = findMaxSpikeNeighbor(first_spike);
    if(Parameters::spikes_to_be_processed.size() != 0) {
        filterAllNeighbors(max_spike, filteredsp);
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

void filterAllNeighbors(Spike max_spike, ofstream& filteredsp) {
    /*Filters all neighbors with smaller amplitudes than max spike that occur
    in the frame of the event.

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
        if(curr_amp < max_spike.amplitude) {
            if(Parameters::verbose) {
                filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << "  " << endl;
            }
            it = Parameters::spikes_to_be_processed.erase(it);
        }
        else {
            //Neighbor has larger amplitude that max neighbor (probably new spike), filter later
            ++it;
        }
    }
    else {
        //Not a neighbor, filter later
        ++it;
    }
	}
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
        curr_frame = curr_spike.frame;
        if(curr_frame <= max_spike.frame + Parameters::noise_duration) {
    		if(areNeighbors(max_spike.channel, curr_channel)) {
                if(!isInnerNeighbor(max_spike.channel, curr_channel)) {
                    if(filteredOuterSpike(curr_spike, max_spike)) {
                        //filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " PN ratio: " << posToNegRatio(curr_spike) << " Area/Amp: " << areaUnderSpike(curr_spike) << " RP time: " << repolarizationTime(curr_spike) << " Filtered by " << max_spike.channel << endl;
                        if(Parameters::verbose) {
                            filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << "  " << endl;
                        }
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
        else {
            //Too far away,filter later
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
    //bool IS_INNER_EVENT = true;
    bool shares_inner = false;
    int NOT_VALID_FRAME = -1;


    int curr_inner_neighbor;
    for(int i = 0; i < Parameters::max_neighbors; i++) {
        curr_inner_neighbor = Parameters::inner_neighbor_matrix[outer_spike.channel][i];
        if(curr_inner_neighbor == -1) {
            break;
        }
        else {
            if(isInnerNeighbor(max_spike.channel, curr_inner_neighbor)) {
                deque<Spike>::iterator it;
                it = Parameters::spikes_to_be_processed.begin();
                //Find max
            	while(it != Parameters::spikes_to_be_processed.end())
            	{
            		if(curr_inner_neighbor ==  it->channel) {
                        if(outer_spike.amplitude < it->amplitude*Parameters::noise_amp_percent) {
                            if(outer_spike.frame < it->frame - Parameters::noise_duration) {
                                //outer spike occurs too far before inner spike, probably new spike
                                shares_inner = true;
                                break;
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
                            break;
                        }
                    }
                    ++it;
            	}
                if(shares_inner == true) {
                    return filtered_spike;
                }
            }
        }
    }
    //int closest_inner_neighbor_channel = getClosestInnerNeighborChannel(outer_spike.channel, max_spike.channel);
    //for(int i = 0; i < Parameters::max_neighbors; i++) {
        float curr_dist;
        float outer_dist_from_center = channelsDist(outer_spike.channel, max_spike.channel);
        for(int i = 0; i < Parameters::max_neighbors; i++) {
            int curr_inner_channel = Parameters::inner_neighbor_matrix[outer_spike.channel][i];
            //out of inner channels
            if(curr_inner_channel == -1) {
                break;
            }
            else {
                curr_dist = channelsDist(curr_inner_channel, max_spike.channel);
                if(curr_dist < outer_dist_from_center) {
                    if(Parameters::masked_channels[curr_inner_channel] == 0) {
                        //masked channel so we relax decay constraints and just return it
                        if(outer_spike.amplitude >= max_spike.amplitude*Parameters::noise_amp_percent) {
                            break;
                        }
                        else {
                            filtered_spike = true;
                            break;
                        }
                    }
                    else {
                        Spike inner_spike = getSpikeFromChannel(curr_inner_channel);
                        if(inner_spike.frame == NOT_VALID_FRAME) {
                            //not a spike, keep searching
                        }
                        else {
                            if(outer_spike.amplitude >= inner_spike.amplitude*Parameters::noise_amp_percent) {
                                //keep searching
                            }
                            else {
                                //Find the closest neighbor in the spike vector and recurse with it (Recursive Step)
                                if(outer_spike.frame <  inner_spike.frame - Parameters::noise_duration) {
                                    //Occured too far before the inner spike to be related
                                    //keep searching
                                }
                                else {
                                    return filteredOuterSpike(inner_spike, max_spike);
                                }
                            }
                        }
                    }
                }
            }
        }
        return filtered_spike;
    //}
}

int getClosestInnerNeighborChannel(int outer_channel, int central_channel) {
    float curr_dist;
    int closest_inner_channel;
    int closest_dist = INT_MAX;
    for(int i = 0; i < Parameters::max_neighbors; i++) {
        int curr_inner_channel = Parameters::inner_neighbor_matrix[outer_channel][i];
        if(curr_inner_channel == -1) {
            break;
        }
        else {
            curr_dist = channelsDist(curr_inner_channel, central_channel);
            if(curr_dist < closest_dist) {
                curr_dist = closest_dist;
                closest_inner_channel = curr_inner_channel;
            }
        }
    }
    return closest_inner_channel;
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
        if(curr_frame <= max_spike.frame + Parameters::noise_duration) {
    		if(areNeighbors(max_spike.channel, curr_channel)) {
                if(isInnerNeighbor(max_spike.channel, curr_channel)) {
                    if(curr_amp < max_spike.amplitude) {
                        if(curr_amp >=  max_spike.amplitude*Parameters::noise_amp_percent) {
                            //filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " PN ratio: " << posToNegRatio(curr_spike)  << " Area/Amp:: " << areaUnderSpike(curr_spike) << " RP time: " << repolarizationTime(curr_spike) << " Filtered by " << max_spike.channel << endl;
                            if(Parameters::verbose) {
                                filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << endl;
                            }
                            it = Parameters::spikes_to_be_processed.erase(it);
                        }
                        else {
                            //filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " PN ratio: " << posToNegRatio(curr_spike) << " Area/Amp:: " << areaUnderSpike(curr_spike) << " RP time: " << repolarizationTime(curr_spike) <<  " Filtered by " << max_spike.channel << endl;
                            if(Parameters::verbose) {
                                filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << endl;
                            }
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
        else {
            //Inner spike occurred much later that max_spike and is only slightly smaller (probably new spike), filter later
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

Spike getSpikeFromChannel(int channel) {
    deque<Spike>::iterator it;
    it = Parameters::spikes_to_be_processed.begin();
    //Find max
    while(it != Parameters::spikes_to_be_processed.end())
    {
        if(channel ==  it->channel) {
            return *it;
        }
        ++it;
    }
    Spike no_spike;
    no_spike.channel = -1;
    no_spike.frame = -1;
    no_spike.amplitude = -1;
    return no_spike;
}


}
