#include "ProcessSpikes.h"

namespace ProcessSpikes {

void filterSpikes(ofstream& spikes_filtered_file,  ofstream& filteredsp)
{
	/*Calls a method from FilterSpikes to filter all spikes. It removes duplicate events
	and writes (or prints) out the channel number, frame, amplitude, and waveforms of the spike.

	Parameters
	----------
	spike_to_be_filtered_file: ofstream&
		The address to the file to be written out to passed in by opened and closed by ProcessSpikes.
	*/
    Spike first_spike = Parameters::spikes_to_be_processed.front();
    Spike max_spike;
	bool isProcessed = false;

	while(!isProcessed) {
        max_spike = FilterSpikes::filterSpikes(first_spike, filteredsp);
		int32_t msc = (int32_t) max_spike.channel;
		int32_t msf = (int32_t) max_spike.frame;
		int32_t msa = (int32_t) max_spike.amplitude;
		int32_t X = (int32_t) 0;
		int32_t Y = (int32_t) 0;

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
				first_spike = Parameters::spikes_to_be_processed.front();
			}
		}
	}
}

void filterLocalizeSpikes(ofstream& spikes_filtered_file, ofstream& filteredsp)
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
    Spike max_spike;
	bool isProcessed = false;

	while(!isProcessed) {

        max_spike = FilterSpikes::filterSpikes(first_spike, filteredsp);

		tuple<float,float> position = LocalizeSpikes::localizeSpike(max_spike);
        if(Parameters::verbose) {
            filteredsp << max_spike.channel << " " << max_spike.frame <<  " " << max_spike.amplitude << "  " << get<0>(position)  << " " << get<1>(position) << endl;
            filteredsp << "E " << Parameters::event_number << endl;
            ++Parameters::event_number;
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

        if(X < 0 || Y < 0) {
            cout << "X real: " << get<0>(position)  << endl;
            cout << "Y real: " << get<1>(position)  << endl;
        }

		if(Parameters::spikes_to_be_processed.size() == 0) {
			isProcessed = true;
		}
		else {
			max_spike = Parameters::spikes_to_be_processed.front();
			if(max_spike.frame > first_spike.frame + Parameters::noise_duration) {
				isProcessed = true;
			}
			else {
				first_spike = max_spike;
			}
		}
	}
}
}
