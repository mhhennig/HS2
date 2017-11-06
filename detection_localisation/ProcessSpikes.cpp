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
	Spike first_spike = Parameters::spikes_to_be_processed.front();
	Parameters::spikes_to_be_processed.pop_front();
	Spike max_spike;
	max_spike = first_spike;
	bool isProcessed = false;

	while(!isProcessed) {
		max_spike = FilterSpikes::filterSpikes(max_spike);
		//cout << "Max Spike: " << max_spike.channel << " " << max_spike.frame << " " << max_spike.amplitude << endl;
		int32_t msc = (int32_t) max_spike.channel;
		int32_t msf = (int32_t) max_spike.frame;
		int32_t msa = (int32_t) max_spike.amplitude;

		spikes_filtered_file.write((char *)&msc, sizeof(msc));
		spikes_filtered_file.write((char *)&msf, sizeof(msf));
		spikes_filtered_file.write((char *)&msa, sizeof(msa));
		spikes_filtered_file.write((char*)&max_spike.written_cutout[0], max_spike.written_cutout.size() * sizeof(int32_t));
		
		if(Parameters::spikes_to_be_processed.size() == 0) {
			isProcessed = true;
		}
		else {
			max_spike = Parameters::spikes_to_be_processed.front();
			if(max_spike.frame > first_spike.frame + Parameters::noise_duration) {
				isProcessed = true;
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

	Spike first_spike = Parameters::spikes_to_be_processed.front();
	Parameters::spikes_to_be_processed.pop_front();
	Spike max_spike;
	max_spike = first_spike;
	bool isProcessed = false;

	while(!isProcessed) {
		max_spike = FilterSpikes::filterSpikes(max_spike);
		//cout << "Max Spike: " << max_spike.channel << " " << max_spike.frame << " " << max_spike.amplitude << endl;
		tuple<float,float> position = LocalizeSpikes::localizeSpike(max_spike);
		//WriteOutSpike to_be_written = {(int32_t)max_spike.channel, (int32_t)max_spike.frame, (int32_t)max_spike.amplitude, get<0>(position), get<1>(position), max_spike.written_cutout};
		//spikes_filtered_file. << max_spike.channel << " " << max_spike.frame << " " << max_spike.amplitude << " " << get<0>(position) << " " << get<1>(position) << " " << cutout.str() <<'\n';
		if(max_spike.frame < 20) {
			stringstream cutout;
			copy(max_spike.written_cutout.begin(), max_spike.written_cutout.end(), ostream_iterator<int>(cutout, " "));
			cout << max_spike.channel << " " << max_spike.frame << " " << max_spike.amplitude << " " << get<0>(position) << " " << get<1>(position) << " " << cutout.str() <<'\n';
		}
		int32_t msc = (int32_t) max_spike.channel;
		int32_t msf = (int32_t) max_spike.frame;
		int32_t msa = (int32_t) max_spike.amplitude;
		int32_t X = (int32_t) floor(get<0>(position) * 1000 + .5);
		int32_t Y = (int32_t) floor(get<1>(position) * 1000 + .5);

		spikes_filtered_file.write((char *)&msc, sizeof(msc));
		spikes_filtered_file.write((char *)&msf, sizeof(msf));
		spikes_filtered_file.write((char *)&msa, sizeof(msa));
		spikes_filtered_file.write((char *)&X, sizeof(X));
		spikes_filtered_file.write((char *)&Y, sizeof(Y));
		spikes_filtered_file.write((char*)&max_spike.written_cutout[0], max_spike.written_cutout.size() * sizeof(int32_t));

		if(Parameters::spikes_to_be_processed.size() == 0) {
			isProcessed = true;
		}
		else {
			max_spike = Parameters::spikes_to_be_processed.front();
			if(max_spike.frame > first_spike.frame + Parameters::noise_duration) {
				isProcessed = true;
			}
			else {
				Parameters::spikes_to_be_processed.pop_front();
			}
		}
	}
}
}