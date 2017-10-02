#include "ProcessSpikes.h"

namespace ProcessSpikes {

void filterSpikes(ofstream& spikes_filtered_file) 
{	
	/*Calls a method from FilterSpikes to filter all spikes. It removes duplicate events
	and writes (or prints) out the channel number, frame, amplitude, and waveforms of the spike.

	Parameters
	----------
	spike_to_be_filtered_file: ofstream&
		The address to the file to be written out to passed in by opened and closed by ProcessSpikes.
	*/
	int curr_channel, curr_frame, curr_amp;
	deque<Spike>::iterator it;
	bool isFinished = false;
	Spike first_spike = Parameters::spikes_to_be_processed.front();
	Parameters::spikes_to_be_processed.pop_front();
	int frame_to_be_filtered = first_spike.frame;
	Spike curr_original_spike;
	curr_original_spike = first_spike;

	while(!isFinished) {
		bool isOriginalSpikeFound = false;
		while(!isOriginalSpikeFound) {
			curr_original_spike = FilterSpikes::filterSpikes(curr_original_spike);
			isOriginalSpikeFound = true;
			deque<Spike>::iterator it;
			it = Parameters::spikes_to_be_processed.begin();
			while(it != Parameters::spikes_to_be_processed.end())
			{
				curr_channel = it->channel;
				curr_amp = it->amplitude;
				for(int i = 0; i < Parameters::max_neighbors; i++) {
					if(Parameters::neighbor_matrix[curr_original_spike.channel][i] == curr_channel) {
						if(curr_amp < curr_original_spike.amplitude) {
							isOriginalSpikeFound = false;
							break;
						}
					}
				}
				++it; 			
			}
		}
		stringstream cutout;
		copy(curr_original_spike.written_cutout.begin(), curr_original_spike.written_cutout.end(), ostream_iterator<int>(cutout, " "));
		spikes_filtered_file << curr_original_spike.channel << " " << curr_original_spike.frame << " " << curr_original_spike.amplitude << " " << cutout.str() << '\n';
		//cout << "Filtred Spike: " << curr_original_spike.channel << " " << curr_original_spike.frame << " " << curr_original_spike.amplitude << '\n';
		Parameters::filtered_spikes += 1;
		if(Parameters::spikes_to_be_processed.size() == 0) {
			isFinished = true;
		}
		else {
			curr_original_spike = Parameters::spikes_to_be_processed.front();
			curr_frame = curr_original_spike.frame;
			if(curr_frame > frame_to_be_filtered + Parameters::noise_duration) {
				isFinished = true;
			}
			else {
				Parameters::spikes_to_be_processed.pop_front();
			}
		}
	}
}

void filterLocalizeSpikes(ofstream& spikes_filtered_file) 
{	
	/*Calls methods from FilterSpikes and LocalizeSpikes to filter and localize
	 all spikes. removing duplicate events and estimating the X and Y coordinates 
	of where they occur. Writes (or prints) out the channel number, frame, amplitude, positions, 
	and waveforms of the spike.

	Parameters
	----------
	spike_to_be_filtered_file: ofstream&
		The address to the file to be written out to passed in by opened and closed by ProcessSpikes.
	*/

	int curr_channel, curr_frame, curr_amp;
	deque<Spike>::iterator it;
	bool isFinished = false;
	Spike first_spike = Parameters::spikes_to_be_processed.front();
	Parameters::spikes_to_be_processed.pop_front();
	int first_frame;
	int frame_to_be_filtered = first_spike.frame;
	Spike curr_original_spike;
	curr_original_spike = first_spike;

	while(!isFinished) {
		bool isOriginalSpikeFound = false;
		first_frame = curr_original_spike.frame;
		while(!isOriginalSpikeFound) {
			//cout << "Removing Dups" << '\n';
			curr_original_spike = FilterSpikes::filterSpikes(curr_original_spike);
			//cout << "Done removing Dups" << '\n';
			isOriginalSpikeFound = true;
			deque<Spike>::iterator it;
			it = Parameters::spikes_to_be_processed.begin();
			while(it != Parameters::spikes_to_be_processed.end())
			{
				curr_channel = it->channel;
				curr_amp = it->amplitude;
				//Checks to see if any duplicate events still exist in the deque
				for(int i = 0; i < Parameters::max_neighbors; i++) {
					if(Parameters::neighbor_matrix[curr_original_spike.channel][i] == curr_channel) {
						if(curr_amp < curr_original_spike.amplitude) {
							isOriginalSpikeFound = false;
							break;
						}
					}
				}
				++it; 			
			}
		}
		tuple<float,float> position = LocalizeSpikes::localizeSpike(curr_original_spike, curr_original_spike.frame - first_frame);
		//Parameters::filtered_spikes += 1;
		//cout << "Spikes Wrten Out: " <<  Parameters::filtered_spikes << "in " << t << " clicks" <<  '\n';
		stringstream cutout;
		copy(curr_original_spike.written_cutout.begin(), curr_original_spike.written_cutout.end(), ostream_iterator<int>(cutout, " "));
		spikes_filtered_file << curr_original_spike.channel << " " << curr_original_spike.frame << " " << curr_original_spike.amplitude << " " << get<0>(position) << " " << get<1>(position) << " " << cutout.str() <<'\n';
		//cout << "Filtered Spike: " << curr_original_spike.channel << " " << curr_original_spike.frame << " " << curr_original_spike.amplitude << '\n';
		if(Parameters::spikes_to_be_processed.size() == 0) {
			isFinished = true;
		}
		else {
			curr_original_spike = Parameters::spikes_to_be_processed.front();
			curr_frame = curr_original_spike.frame;
			if(curr_frame > frame_to_be_filtered + Parameters::noise_duration) {
				isFinished = true;
			}
			else {
				Parameters::spikes_to_be_processed.pop_front();
			}
		}
	}
}
}