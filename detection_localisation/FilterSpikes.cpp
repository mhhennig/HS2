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
    cout << "Filtering: " << first_spike.channel << " " << first_spike.frame <<  " " << first_spike.amplitude << endl;
    Spike max_spike;
    cout << "Finding Max Spike..." << endl;
    max_spike = findMaxSpikeNeighbor(first_spike);
    cout << "Max Spike: " << max_spike.channel << " " << max_spike.frame <<  " " << max_spike.amplitude << endl;
    if(Parameters::spikes_to_be_processed.size() != 0) {
        cout << "Updating neighbors max spike..." << endl;
        max_spike = updateNeighborsMaxSpike(max_spike);
        cout << "Updated!" << endl;
        cout << "Filtering outer neighbors..." << endl;
        //Find the max amplitude neighbor of the first spike
        filterOuterNeighbors(max_spike, filteredsp);
        cout << "Filtered outer neighors!" << endl;
        cout << "Filtering inner neighbors..." << endl;
        filterInnerNeighbors(max_spike, filteredsp);
        cout << "Filtered inner neighors!" << endl;
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
    cout << index << endl;
    cout << Parameters::spikes_to_be_processed.size() << endl;
    cout << max_spike.channel << endl;
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
		The max spike detected in the event.

	*/

    int curr_channel, curr_amp;
    Spike curr_spike;
    vector<tuple<int, int>> inner_spiking_neighbors;
    vector<tuple<int, int>> outer_spiking_neighbors;

    //Find the max amplitude neighbor of the first spike
    deque<Spike>::iterator it;
	it = Parameters::spikes_to_be_processed.begin();
	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		curr_amp = it->amplitude;
        if(areNeighbors(max_spike.channel, curr_channel)) {
            if(channelsDist(max_spike.channel, curr_channel) <= Parameters::inner_radius) {
                inner_spiking_neighbors.push_back(make_tuple(curr_channel, curr_amp));
            }
            else {
                outer_spiking_neighbors.push_back(make_tuple(curr_channel, curr_amp));
            }
	    }
        ++it;
	}

    max_spike = updateOuterNeighbors(max_spike, outer_spiking_neighbors);
    max_spike = updateInnerNeighbors(max_spike, inner_spiking_neighbors);

    return max_spike;
}

Spike updateOuterNeighbors(Spike max_spike, vector<tuple<int, int>> outer_spiking_neighbors) {
    /*Updates all the outer neighbors that spiked

	Parameters
	----------
	max_spike: Spike
        The largest amplitude spike belonging to the event of the first spike.
    */
    vector<tuple<int, int>>::iterator it;
    it = max_spike.outer_neighbors.begin();
    tuple<int, int> curr_neighbor;
    int curr_channel;
    int index = it - max_spike.outer_neighbors.begin();
    while(it != max_spike.outer_neighbors.end())
    {
        curr_neighbor = *it;
        curr_channel = get<0>(curr_neighbor);
        index = it - max_spike.outer_neighbors.begin();

        vector<tuple<int, int>>::iterator it2;
        it2 = outer_spiking_neighbors.begin();
        tuple<int, int> curr_spike_neighbor;
        int curr_spike_channel, curr_spike_amp;
        while(it2 != outer_spiking_neighbors.end())
        {
            curr_spike_neighbor = *it2;
            curr_spike_channel = get<0>(curr_spike_neighbor);
            curr_spike_amp = get<1>(curr_spike_neighbor);
            if(curr_spike_channel == curr_channel) {
                max_spike.outer_neighbors.at(index) = make_tuple(curr_spike_channel, curr_spike_amp);
                break;
            }
            else {
                ++it2;
            }
        }
        ++it;
    }
    return max_spike;

}

Spike updateInnerNeighbors(Spike max_spike, vector<tuple<int, int>> inner_spiking_neighbors) {
    /*Updates all the inner neighbors that spiked

	Parameters
	----------
	max_spike: Spike
        The largest amplitude spike belonging to the event of the first spike.
    */

    vector<tuple<int, int>>::iterator it;
	it = max_spike.inner_neighbors.begin();
    tuple<int, int> curr_neighbor;
    int curr_channel;
    int index = 0;
	while(it != max_spike.inner_neighbors.end())
	{
		curr_neighbor = *it;
        curr_channel = get<0>(curr_neighbor);
        index = it - max_spike.inner_neighbors.begin();

        vector<tuple<int, int>>::iterator it2;
    	it2 = inner_spiking_neighbors.begin();
        tuple<int, int> curr_spike_neighbor;
        int curr_spike_channel, curr_spike_amp;
    	while(it2 != inner_spiking_neighbors.end())
    	{
    		curr_spike_neighbor = *it2;
            curr_spike_channel = get<0>(curr_spike_neighbor);
            curr_spike_amp = get<1>(curr_spike_neighbor);
            if(curr_spike_channel == curr_channel) {
                max_spike.inner_neighbors.at(index) = make_tuple(curr_spike_channel, curr_spike_amp);
                break;
            }
            else {
                ++it2;
            }
    	}
        ++it;
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
	int curr_channel;
    deque<Spike>::iterator it;
	it = Parameters::spikes_to_be_processed.begin();

	while(it != Parameters::spikes_to_be_processed.end())
	{
		curr_spike = *it;
		curr_channel = it->channel;
		if(areNeighbors(max_spike.channel, curr_channel)) {
            if(channelsDist(max_spike.channel, curr_channel) > Parameters::inner_radius) {
                if(filteredOuterSpike(curr_spike, max_spike)) {
                    filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " " << "Filtered by " << max_spike.channel << endl;
                    it = Parameters::spikes_to_be_processed.erase(it);
                }
                else {
                    //Not a decaying outer spike (probably a new spike), filter later
                    ++it;
                }
            }
            else {
                //Inner neighbor - all inner neighbors that needed filtering should have been filtered earlier, filter later
                ++it;
            }
		}
		else {
            //Not a neighbor, filter later
			++it;
		}
	}
}

bool filteredOuterSpike(Spike outer_spike, Spike max_spike) {
    //True if should be filtered, false by default
    bool filtered_spike = false;

    vector<tuple<int, int>>::iterator it;
    it = outer_spike.inner_neighbors.begin();
    //Searching for shared inner neighbor (Base Case search)
    tuple<int, int> curr_inner_neighbor;
    int curr_inner_channel, curr_inner_amp;
    while(it != outer_spike.inner_neighbors.end())
    {
        curr_inner_neighbor = *it;
        curr_inner_channel = get<0>(curr_inner_neighbor);
        curr_inner_amp = get<1>(curr_inner_neighbor);

        vector<tuple<int, int>>::iterator it2;
        it2 = max_spike.inner_neighbors.begin();
        tuple<int, int> curr_max_inner_neighbor;
        int curr_max_inner_channel, curr_max_inner_amp;
        while(it2 != max_spike.inner_neighbors.end())
        {
            curr_max_inner_neighbor = *it2;
            curr_max_inner_channel = get<0>(curr_max_inner_neighbor);
            curr_max_inner_amp = get<1>(curr_max_inner_neighbor);
            //Shares an inner channel with the max spike (Base Case)
            if(curr_inner_channel == curr_max_inner_channel) {
                if(outer_spike.amplitude < curr_max_inner_amp*Parameters::noise_amp_percent) {
                    cout <<"Outer spike filter: " << outer_spike.channel << endl;
                    cout << curr_max_inner_channel << endl;
                    cout << curr_max_inner_amp*Parameters::noise_amp_percent << endl;
                    int inner_spike_frame = 0;
                    int curr_channel, curr_frame;
                    Spike curr_spike;
                    deque<Spike>::iterator it3;
                    it3 = Parameters::spikes_to_be_processed.begin();
                    //Get frame of shared inner neighbor spike
                    while(it3 != Parameters::spikes_to_be_processed.end())
                    {
                        curr_spike = *it3;
                        curr_channel = it3->channel;
                        curr_frame = it3->frame;
                        if(curr_channel == curr_inner_channel) {
                            inner_spike_frame = curr_frame;
                            break;
                        }
                        else {
                            ++it3;
                        }
                    }
                    if(outer_spike.frame < inner_spike_frame - Parameters::noise_duration) {
                        //outer spike occurs too far before inner spike, probably new spike
                        return filtered_spike;
                    }
                    else {
                        filtered_spike = true;
                        return filtered_spike;
                    }
                }
                else {
                    //amplitude too big to be a duplicate spike
                    return filtered_spike;
                }
            }
            else {
                ++it2;
            }
        }
        ++it;
    }
    //Doesn't share an inner_neighbor with max_spike (find closest inner neighbor)
    //INT_MAX
    int closest_dist = 10000;
    int closest_channel = -1;
    int closest_amplitude = 0;
    int curr_dist;
    it = outer_spike.inner_neighbors.begin();
    while(it != outer_spike.inner_neighbors.end())
    {
        curr_inner_neighbor = *it;
        curr_inner_channel = get<0>(curr_inner_neighbor);
        curr_inner_amp = get<1>(curr_inner_neighbor);
        curr_dist = channelsDist(curr_inner_channel, max_spike.channel);
        if(curr_dist < closest_dist) {
            closest_dist = curr_dist;
            closest_channel = curr_inner_channel;
            closest_amplitude = curr_inner_amp;
        }
    }

    if(outer_spike.amplitude >= closest_amplitude*Parameters::noise_amp_percent) {
        return filtered_spike;
    }
    else {
        //Find the closest neighbor in the spike vector and recurse with it (Recursive Step)
        int curr_channel, curr_frame;
        Spike curr_spike;
        deque<Spike>::iterator it;
    	it = Parameters::spikes_to_be_processed.begin();

    	while(it != Parameters::spikes_to_be_processed.end())
    	{
    		curr_spike = *it;
    		curr_channel = it->channel;
            curr_frame = it->frame;
    		if(curr_channel == closest_channel) {
                if(outer_spike.frame >=  curr_frame - Parameters::noise_duration) {
                    return filteredOuterSpike(curr_spike, max_spike);
                }
                else {
                    //Occured too far before the inner spike to be related
                    return filtered_spike;
                }
    		}
    		else {
    			++it;
    		}
    	}
    }
    return filtered_spike;
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
            if(channelsDist(max_spike.channel, curr_channel) <= Parameters::inner_radius) {
                if(curr_amp < max_spike.amplitude) {
                    if(curr_amp >=  max_spike.amplitude*Parameters::noise_amp_percent) {
                        if(curr_frame >  max_spike.frame + Parameters::noise_duration) {
                            //Inner spike occurred much later that max_spike and is only slightly smaller (probably new spike), filter later
                            ++it;
                        }
                        else {
                            filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " " << "Filtered by " << max_spike.channel << endl;
                            it = Parameters::spikes_to_be_processed.erase(it);
                        }
                    }
                    else {
                        filteredsp << curr_spike.channel << " " << curr_spike.frame <<  " " << curr_spike.amplitude << " " << "Filtered by " << max_spike.channel << endl;
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

}
