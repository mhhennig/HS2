#include "SpikeHandler.h"

int Parameters::num_channels;
int Parameters::num_recording_channels;
int Parameters::spike_delay;
int Parameters::spike_peak_duration;
int Parameters::noise_duration;
float Parameters::noise_amp_percent;
int Parameters::max_neighbors;
int** Parameters::neighbor_matrix;
int** Parameters::channel_positions;
int Parameters::aGlobal;
bool Parameters::to_localize;
int** Parameters::baselines;
int Parameters::cutout_start;
int Parameters::cutout_end;
int Parameters::index_data;
int Parameters::index_baselines;
int Parameters::iterations;
int Parameters::frames;
int Parameters::maxsl;
int Parameters::end_raw_data;
int Parameters::bad_index = 0;
short* Parameters::raw_data;
deque<Spike> Parameters::spikes_to_be_processed;
std::ofstream spikes_filtered_file;


void setInitialParameters(int _num_channels, int _num_recording_channels, int _spike_delay, int _spike_peak_duration, string file_name, \
						  int _noise_duration, float _noise_amp_percent, int** _channel_positions, int** _neighbor_matrix, \
						  int _max_neighbors, bool _to_localize = false, int _cutout_start= 10, int _cutout_end=20, int _maxsl = 0) 
{
	/*This sets all the initial parameters needed to run the filtering algorithm.

	Parameters
	----------
	_num_channels: int
		Number of channels on the probe.
	_num_recording_channels: int
		Number of channels to be used for spike data.
	_spike_delay: int
		The number of frames back a spike occurred after it was detected (where the beginning of the spike was).
	_spike_peak_duration: int
		The number of frames it takes a spike amplitude to fully decay.
	_noise_duration: int
		The frames in which the true spike can occur after the first detection of the spike. 
		Ideally, the _noise_duration would be zero if there is no noise (then the first spike 
		that occurs is the original spike), but sometimes the true spike is detected after
		a duplicate.
	_noise_amp_percent: float
		The amplitude percent difference at which two spikes can be considered unique 
		even if the spike detected after has a smaller amplitude (with zero _noise_amp_percent, 
		two concurrent spikes where the second spike has a slightly smaller amplitude 
		will be considered duplicates).
	_channel_positions: 2D int array
		Indexed by the channel number starting at 0 and going up to num_recording_channels - 1. Each 
		index contains pointer to another array which contains X and Y position of the channel. User creates 
		this before calling SpikeHandler. 
	_neighbor_matrix: 2D int array
		Indexed by the channel number starting at 0 and going up to num_recording_channels - 1. Each 
		index contains pointer to another array which contains channel number of all its neighbors.
		User creates this before calling SpikeHandler.
	_max_neighbors: int
		The maximum number of neighbors a channel can have (a neighbor is any channel that can receive the
		same spike waveform as the original channel if the spike occurred at the original channel).
	_to_localize: bool
		True: Localize the spike using our localization method
		False: Do not localize.
	cutout_length: int
		The cutout length to be written out for the spike. Can't be larger than extra data tacked on to raw data.
	_spikes_to_be_processed: Spike deque
		Contains all spikes to be proccessed when spike_peak_duration number of frames is stored.
	_maxsl: int
		The number of frames until a spike is accepted when the peak value is given.
	*/
	if(_num_channels < 0 || _num_channels < _num_recording_channels) {
		cout << "Number of channels given incorrectly. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_num_recording_channels < 0) {
		cout << "Number of recording channels less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_max_neighbors < 0) {
		cout << "Number of max neighbors less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_spike_delay < 0) {
		cout << "Spike Delay less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_spike_peak_duration < 0) {
		cout << "Spike Peak Duration less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_noise_amp_percent < 0 || _noise_amp_percent > 1) {
		cout << "Noise Amplitude Percent not a valid percentage. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_noise_duration < 0) {
		cout << "Cutout Length less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_cutout_start < 0) {
		cout << "Cutout Start less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_cutout_end < 0) {
		cout << "Cutout End less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_maxsl < 0) {
		cout << "Maxsl less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	Parameters::num_channels = _num_channels;
	Parameters::num_recording_channels = _num_recording_channels;
	Parameters::max_neighbors = _max_neighbors;
	Parameters::spike_delay = _spike_delay;
	Parameters::spike_peak_duration = _spike_peak_duration;
	Parameters::noise_duration = _noise_duration;
	Parameters::noise_amp_percent = _noise_amp_percent;
	Parameters::to_localize = _to_localize;
	Parameters::channel_positions = _channel_positions;
	Parameters::neighbor_matrix = _neighbor_matrix;
	Parameters::cutout_start = _cutout_start;
	Parameters::cutout_end = _cutout_end;
	Parameters::maxsl = _maxsl;

	spikes_filtered_file.open(file_name);


}
void loadRawData(short* _raw_data, int _index_data, int _iterations, int _frames, int _additional_data) {
	/*Every iteration where new raw data is passed in, this sets pointer to new data and gives the
	index to start accessing the data at

	Parameters
	----------
	_raw_data: short array
		Stores all the raw data for the iteration. The data is formatted where the current reading for the first frame
		at every channel (including non-spike detecting channels) is given in order from smallest channel to largest
		by channel number. Then it gives the reading for the second frame at every channel. For each iteration, buffer
		data should be provided for the cutout.
	_index_data: int
		The index at which the readings for the data start. This allows the spike detection to skip the buffer data used
		only for cutout data.
	_iterations: int
		Number of current iterations of raw data passed in. User starts this at 0 and increments it for each chunk of data
		passed into loadRawData.
	_frames: int
		The number of frames passed into loadRawData EXCLUDING the buffer frames.
	*/
	if(_index_data < 0) {
		spikes_filtered_file.close();
		cout << "Index Data less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_iterations < 0) {
		spikes_filtered_file.close();
		cout << "Iterations less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_frames < 0) {
		spikes_filtered_file.close();
		cout << "Input frames less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	if(_additional_data < 0) {
		spikes_filtered_file.close();
		cout << "Additional data less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}

	Parameters::raw_data = _raw_data;
	Parameters::index_data = _index_data;
	Parameters::iterations = _iterations;
	Parameters::frames = _frames;
	Parameters::end_raw_data = (Parameters::frames + _additional_data + Parameters::index_data)*Parameters::num_channels + Parameters::num_channels - 1;
}
void setLocalizationParameters(int _aGlobal, int** _baselines, int _index_baselines) {
	/*Sets all time dependent variables for localization. The localization needs a global noise value
	and a precalculated median baseline value for processing amplitudes of spikes.

	Parameters
	----------
	_aGlobal: int
		The global noise values for all the channels.
	_baselines: int
		Contains spike_delay number of frames of median baseline values. Updated by user at every frame.
	*/
	if(_index_baselines < 0) {
		spikes_filtered_file.close();
		cout << "Index baselines less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}
	Parameters::aGlobal = _aGlobal;
	Parameters::baselines = _baselines;
	Parameters::index_baselines = _index_baselines;
}

void addSpike(int channel, int frame, int amplitude) {
	/*Adds a spike to the spikes_to_be_processed deque. Once the frame of the spike to be added is
	greater than the spike_peak_duration larger than the first spike in the deque, it will process 
	all the current spikes in the deque and then attempt to add the spike again. It will keep repeating 
	this process until the spike to be added frame is smaller than the spike_peak_duration larger than the
	first spike or the deque is empty.

	Parameters
	----------
	channel: int
		The channel at which the spike is detected.
	frame: int
		The frame at which the spike is detected.
	amplitude: int
		The amplitude at which the spike is detected.
	*/
	int cutout_size = Parameters::cutout_start + Parameters::cutout_end + 1;
	int amp_cutout_size = Parameters::spike_delay*2 + 1;
	int frames_processed = Parameters::frames*Parameters::iterations;

	if(channel < Parameters::num_recording_channels && channel >= 0) {
		int curr_neighbor_channel, curr_reading;
		Spike spike_to_be_added;
		spike_to_be_added.channel = channel;
		spike_to_be_added.frame = frame;
		spike_to_be_added.amplitude = amplitude;
		int ASCALE = -64;

		for(int i = 0; i < cutout_size; i++) {
			try {
				int curr_reading_index = (frame - Parameters::cutout_start - frames_processed + Parameters::index_data + i)*Parameters::num_channels + channel;
				if(curr_reading_index < 0 || curr_reading_index > Parameters::end_raw_data) {
					Parameters::bad_index += 1;
					curr_reading = 0;
				}
				else {
					curr_reading = Parameters::raw_data[curr_reading_index];
				}
			} catch (...) {
				spikes_filtered_file.close();
				cout << "Raw Data and it parameters entered incorrectly, could not access data. Terminating SpikeHandler." << endl;
				exit(EXIT_FAILURE);
			}
			spike_to_be_added.written_cutout.push_back(curr_reading);
}

		for(int i = 0; i < Parameters::max_neighbors; i++) {
			try {
				curr_neighbor_channel = Parameters::neighbor_matrix[channel][i];
			} catch (...) { 
				spikes_filtered_file.close();
				cout << "Neighbor matrix improperly created. Terminating SpikeHandler" << endl;
				exit(EXIT_FAILURE);
			} 
			if(curr_neighbor_channel != -1) {
				for(int j = 0; j < amp_cutout_size; j++) {
					try {
  						curr_reading = Parameters::raw_data[(frame - Parameters::spike_delay - frames_processed + Parameters::index_data + j)*Parameters::num_channels + curr_neighbor_channel];
					} catch (...) { 
						spikes_filtered_file.close();
						cout << "Raw Data and it parameters entered incorrectly, could not access data. Terminating SpikeHandler." << endl;
						exit(EXIT_FAILURE);
					}

					int curr_amp = ((curr_reading - Parameters::aGlobal) * ASCALE - Parameters::baselines[curr_neighbor_channel][Parameters::index_baselines]);
					spike_to_be_added.amp_cutouts.push_back(curr_amp);
				} 
			}
			//Out of neighbors to add cutout for
			else {
				break;
			}
		}
		//Add Spike to the current frame window
		//cout << "Spike to be added: " << spike_to_be_added.channel << " " << spike_to_be_added.frame << " " << spike_to_be_added.amplitude << endl;
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
				if(current_frame > first_frame + (Parameters::spike_peak_duration + Parameters::noise_duration)) {
					if(Parameters::to_localize) {
						try {
							ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file);
							//spikes_filtered_file.close();
							//cout << "doneg" << endl;
							//exit(EXIT_FAILURE);
						} catch (...) { 
							spikes_filtered_file.close();
							cout << "Baseline matrix or its parameters entered incorrectly. Terminating SpikeHandler." << endl;
							exit(EXIT_FAILURE);
						} 
					}
					else {
						ProcessSpikes::filterSpikes(spikes_filtered_file);
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
	//Filter any remaining spikes leftover at the end and close the spike file.
	while(Parameters::spikes_to_be_processed.size() != 0){
		if(Parameters::to_localize) {
			ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file); 
		}
		else {
			ProcessSpikes::filterSpikes(spikes_filtered_file); 
		}	
	}
	cout << "Cutouts that had to be changed to 0 because not in data: " <<  Parameters::bad_index << endl;
	spikes_filtered_file.close();
}