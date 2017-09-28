#include "SpikeHandler.h"

int Parameters::num_channels;
int Parameters::num_recording_channels;
int Parameters::spike_delay;
int Parameters::spike_peak_duration;
int Parameters::noise_duration;
int Parameters::noise_amp;
int Parameters::max_neighbors;
int** Parameters::neighbor_matrix;
int** Parameters::channel_positions;
int Parameters::aGlobal;
bool Parameters::to_localize;
int* Parameters::baselines;
int Parameters::start_cutout;
int Parameters::end_cutout;
int Parameters::filtered_spikes;
deque<int> Parameters::amps;
short* Parameters::raw_data;
deque<Spike> Parameters::spikes_to_be_processed;
std::ofstream spikes_filtered_file;

void setInitialParameters(int _num_channels, int _num_recording_channels, int _spike_delay, int _spike_peak_duration, \
						  int _noise_duration, int _noise_amp, int** _channel_positions, int** _neighbor_matrix, \
						  int _max_neighbors, bool _to_localize, int _start_cutout, int _end_coutout) 
{
	/*This sets all the initial parameters needed to run the filtering algorithm.

	Parameters
	----------
	_num_channels: int
		The number of channels in the probe.
	_num_recording_channels: int
		The number of channels in the probe that are actually recording. Not ground/reference channels.
	_spike_delay: int
		The frames before the spike is detected that the spike actually start. Essential for localization.
	_spike_peak_duration: int
		The frames that a spike amplitude may last before disappearing. This helps find duplicate spikes 
		since they must occur before this duration is over since the first spike fired.
	_noise_duration: int
		The frames in which the true spike can occur after the first detection of the spike. 
		Ideally, the _noise_duration would be zero if there is no noise (then the first spike 
		that occurs is the original spike), but sometimes the true spike is detected after
		a duplicate.
	_noise_amp: int
		The amplitude difference at which two spikes can be considered the duplicates 
		even if the spike detected after has a larger amplitude (with zero _noise_amp, 
		two concurrent spikes where the second spike has a slightly larger amplitude 
		will be considered unique).
	_channel_positions: 2D int array
		A 2D array with each index representing channel numbers containing arrays of 
		size two that hold the X and Y coordinates of the channel. The length of 
		Rows = _num_recording_channels and Cols = 2.
	_neighbor_matrix: 2D int array
		A 2D array with each index representing channel numbers that  correspond to integer 
		arrays that contain the channel numbers that the index channel is neighboring.
		The length of Rows = _num_recording_channels and Cols = _max_neighbors.
	_max_neighbors: int
		The maximum number of neighbors a channel can have (a neighbor is any channel that can receive the
		same spike waveform as the original channel if the spike occurred at the original channel).
	_to_localize: bool
		True: Localize and spike using our localization method
		False: Do not localize.
	_start_cutout: int
		The number of frames before the spike is detected. Used for writing out waveform.
	_end_cutout: int
		The number of frames after the spike is detected. Used for writing out waveform.
	_spikes_to_be_processed: Spike deque
		The current spikes to be processed when a chunk is full. A Spike is a struct defined
		by an channel, frame, and amplitude.

	*/
	Parameters::num_channels = _num_channels;
	Parameters::num_recording_channels = _num_recording_channels;
	Parameters::max_neighbors = _max_neighbors;
	Parameters::spike_delay = _spike_delay;
	Parameters::spike_peak_duration = _spike_peak_duration;
	Parameters::noise_duration = _noise_duration;
	Parameters::noise_amp = _noise_amp;
	Parameters::to_localize = _to_localize;
	Parameters::channel_positions = _channel_positions;
	Parameters::neighbor_matrix = _neighbor_matrix;
	Parameters::start_cutout = _start_cutout;
	Parameters::end_cutout = _end_coutout;
	Parameters::filtered_spikes = 0;
	spikes_filtered_file.open("SpikesFiltered");

}
void loadRawData(short* _raw_data) {
	Parameters::raw_data = _raw_data;
}
void setLocalizationParameters(int _aGlobal, int* _baselines) {
	Parameters::aGlobal = _aGlobal;
	Parameters::baselines = _baselines;
}

void addSpike(int channel, int frame, int amplitude) {
	if(channel < Parameters::num_recording_channels && channel > 0) {
		Spike spike_to_be_added;
		spike_to_be_added.channel = channel;
		spike_to_be_added.frame = frame;
		spike_to_be_added.amplitude = amplitude;
		//cout << "Adding Spike" << '\n';
		bool isAdded = false;
		while(!isAdded) {
			if(Parameters::spikes_to_be_processed.empty()) {
				Parameters::spikes_to_be_processed.push_back(spike_to_be_added);
				isAdded = true;
			}
			else {
				Spike first_spike = Parameters::spikes_to_be_processed.front();
				int first_frame = first_spike.frame;
				int current_frame = spike_to_be_added.frame;
				//cout << spikes_to_be_processed.size() << '\n';
				if(current_frame > first_frame + (Parameters::spike_peak_duration + Parameters::noise_duration)) {
					if(Parameters::to_localize) {
						ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file);
						//cout << "Done filtering Spike" << '\n';
					}
					else {
						//cout << "Filtering Spik" << '\n';
						ProcessSpikes::filterSpikes(spikes_filtered_file);
						//cout << "Done filtering Spike" << '\n';

					}
				} 
				else {
					Parameters::spikes_to_be_processed.push_back(spike_to_be_added);
					isAdded = true;
				}
			}
		}
	}
}
void terminateSpikeHandler() {
	while(Parameters::spikes_to_be_processed.size() != 0){
		cout << "Tryin" << '\n';
		if(Parameters::to_localize) {
			ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file); 
		}
		else {
			ProcessSpikes::filterSpikes(spikes_filtered_file); 
		}	
	}
	spikes_filtered_file << "done" << '\n';
	spikes_filtered_file.close();
}