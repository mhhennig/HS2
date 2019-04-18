#include "SpikeHandler.h"

int Parameters::ASCALE = -64;
int Parameters::num_com_centers;
int Parameters::num_channels;
int Parameters::spike_peak_duration;
int Parameters::noise_duration;
float Parameters::noise_amp_percent;
int Parameters::max_neighbors;
int **Parameters::neighbor_matrix;
float **Parameters::channel_positions;
int **Parameters::inner_neighbor_matrix;
int **Parameters::outer_neighbor_matrix;
int Parameters::aGlobal;
bool Parameters::to_localize;
bool Parameters::verbose;
int **Parameters::baselines;
int Parameters::cutout_start;
int Parameters::cutout_end;
int Parameters::index_data;
int Parameters::index_baselines;
int Parameters::iterations;
int Parameters::max_frames_processed;
int Parameters::before_chunk;
int Parameters::after_chunk;
int Parameters::maxsl;
int Parameters::end_raw_data;
short *Parameters::raw_data;
int *Parameters::masked_channels;
int Parameters::event_number;
bool Parameters::debug;
bool Parameters::decay_filtering;
float Parameters::inner_radius;

deque<Spike> Parameters::spikes_to_be_processed;
std::ofstream filteredsp;
std::ofstream spikes_filtered_file;

namespace SpikeHandler {

void setInitialParameters(int _num_channels,
                          int _spike_peak_duration, string file_name,
                          int _noise_duration, float _noise_amp_percent,
                          float _inner_radius, int *_masked_channels,
                          float **_channel_positions, int **_neighbor_matrix,
                          int _max_neighbors, int _num_com_centers = 1,
                          bool _to_localize = false, int _cutout_start = 10,
                          int _cutout_end = 20, int _maxsl = 0,
                          bool _decay_filtering = true, bool _verbose = true) {
  /*This sets all the initial parameters needed to run the filtering algorithm.

  Parameters
  ----------
  _num_channels: int
          Number of channels on the probe
(where the beginning of the spike was).
  file_name: string
  The name of the file that the processed spikes will be written to in binary.
  _noise_duration: int
          The frames in which the true spike can occur after the first detection
of the spike.
          Ideally, the _noise_duration would be zero if there is no noise (then
the first spike
          that occurs is the original spike), but sometimes the true spike is
detected after
          a duplicate.
  _noise_amp_percent: float
          The amplitude percent difference at which two spikes can be considered
unique
          even if the spike detected after has a smaller amplitude (with zero
_noise_amp_percent,
          two concurrent spikes where the second spike has a slightly smaller
amplitude
          will be considered duplicates).
_masked_channels: 1D int array
  The channels that are masked for the detection (AKA we get input from them)
  _channel_positions: 2D int array
          Indexed by the channel number starting at 0 and going up to
num_channels: Each
          index contains pointer to another array which contains X and Y
position of the channel. User creates
          this before calling SpikeHandler.
  _neighbor_matrix: 2D int array
          Indexed by the channel number starting at 0 and going up to
num_channels: Each
          index contains pointer to another array which contains channel number
of all its neighbors.
          User creates this before calling SpikeHandler.
  _max_neighbors: int
          The maximum number of neighbors a channel can have (a neighbor is any
channel that can receive the
          same spike waveform as the original channel if the spike occurred at
the original channel).
  _to_localize: bool
          True: Localize the spike using our localization method
          False: Do not localize.
  cutout_length: int
          The cutout length to be written out for the spike. Can't be larger
than extra data tacked on to raw data.
  _spikes_to_be_processed: Spike deque
          Contains all spikes to be proccessed when spike_peak_duration number
of frames is stored.
  _maxsl: int
          The number of frames until a spike is accepted when the peak value is
given.
  */
  if (_num_channels < 0) {
    cout << "Number of channels given incorrectly. Terminating Spike Handler"
         << endl;
    exit(EXIT_FAILURE);
  }
  if (_max_neighbors < 0) {
    cout << "Number of max neighbors less than 0. Terminating Spike Handler"
         << endl;
    exit(EXIT_FAILURE);
  }
  if (_spike_peak_duration < 0) {
    cout << "Spike Peak Duration less than 0. Terminating Spike Handler"
         << endl;
    exit(EXIT_FAILURE);
  }
  if (_noise_amp_percent < 0 || _noise_amp_percent > 1) {
    cout << "Noise Amplitude Percent not a valid percentage. Terminating Spike "
            "Handler"
         << endl;
    exit(EXIT_FAILURE);
  }
  if (_noise_duration < 0) {
    cout << "Cutout Length less than 0. Terminating Spike Handler" << endl;
    exit(EXIT_FAILURE);
  }
  if (_cutout_start < 0) {
    cout << "Cutout Start less than 0. Terminating Spike Handler" << endl;
    exit(EXIT_FAILURE);
  }
  if (_cutout_end < 0) {
    cout << "Cutout End less than 0. Terminating Spike Handler" << endl;
    exit(EXIT_FAILURE);
  }
  if (_maxsl < 0) {
    cout << "Maxsl less than 0. Terminating Spike Handler" << endl;
    exit(EXIT_FAILURE);
  }
  if (_inner_radius < 0) {
    cout << "Inner Radius less than 0. Terminating Spike Handler" << endl;
    exit(EXIT_FAILURE);
  }

  Parameters::num_channels = _num_channels;
  Parameters::max_neighbors = _max_neighbors;
  Parameters::spike_peak_duration = _spike_peak_duration;
  Parameters::noise_duration = _noise_duration;
  Parameters::noise_amp_percent = _noise_amp_percent;
  Parameters::to_localize = _to_localize;
  Parameters::channel_positions = _channel_positions;
  Parameters::neighbor_matrix = _neighbor_matrix;
  Parameters::cutout_start = _cutout_start;
  Parameters::cutout_end = _cutout_end;
  Parameters::maxsl = _maxsl;
  Parameters::masked_channels = _masked_channels;
  Parameters::inner_radius = _inner_radius;
  Parameters::event_number = 0;
  Parameters::debug = false;
  Parameters::verbose = _verbose;
  Parameters::decay_filtering = _decay_filtering;
  Parameters::num_com_centers = _num_com_centers;

  Parameters::inner_neighbor_matrix = createInnerNeighborMatrix();
  Parameters::outer_neighbor_matrix = createOuterNeighborMatrix();
  fillNeighborLayerMatrices();
  if (Parameters::debug) {
    for (int i = 0; i < Parameters::num_channels; i++) {
      cout << "Channel: " << i << endl;
      cout << "Inner Neighbors: ";
      for (int j = 0; j < Parameters::max_neighbors; j++) {
        cout << Parameters::inner_neighbor_matrix[i][j] << "  ";
      }
      cout << endl;
      cout << "Outer Neighbors: ";
      for (int k = 0; k < Parameters::max_neighbors; k++) {
        cout << Parameters::outer_neighbor_matrix[i][k] << "  ";
      }
      cout << endl;
    }
  }

  if(spikes_filtered_file.is_open()){
    spikes_filtered_file.close();
  }
//   spikes_filtered_file.open(file_name + ".bin", ios::trunc | ios::binary);
  spikes_filtered_file.open(file_name + ".bin", ios::binary);
  if (_verbose) {
    filteredsp.open(file_name + "_filtered_spikes.asc");
  }
    Parameters::spikes_to_be_processed.clear();
}
void loadRawData(short *_raw_data, int _index_data, int _iterations,
                 int maxFramesProcessed, int before_chunk, int after_chunk) {
  /*Every iteration where new raw data is passed in, this sets pointer to new
  data and gives the
  index to start accessing the data at

  Parameters
  ----------
  _raw_data: short array
          Stores all the raw data for the iteration. The data is formatted where
  the current reading for the first frame
          at every channel (including non-spike detecting channels) is given in
  order from smallest channel to largest
          by channel number. Then it gives the reading for the second frame at
  every channel. For each iteration, buffer
          data should be provided for the cutout.
  _index_data: int
          The index at which the readings for the data start. This allows the
  spike detection to skip the buffer data used
          only for cutout data.
  _iterations: int
          Number of current iterations of raw data passed in. User starts this
  at 0 and increments it for each chunk of data
          passed into loadRawData.
  maxFramesProcessed: int
          The max number of frames passed sorting algorithm
  frames.
  before_chunk: int
          The number of buffer frames before tInc
  after_chunk: int
          The number of buffer frames after tInc
  */
  Parameters::raw_data = _raw_data;
  Parameters::index_data = _index_data;
  Parameters::iterations = _iterations;
  Parameters::max_frames_processed = maxFramesProcessed;
  Parameters::before_chunk = before_chunk;
  Parameters::after_chunk = after_chunk;
  Parameters::end_raw_data =
      (Parameters::max_frames_processed + Parameters::after_chunk  + Parameters::index_data) *
          Parameters::num_channels +
      Parameters::num_channels - 1;
}
void setLocalizationParameters(int _aGlobal, int **_baselines,
                               int _index_baselines) {
  /*Sets all time dependent variables for localization. The localization needs a
  global noise value
  and a precalculated median baseline value for processing amplitudes of spikes.

  Parameters
  ----------
  _aGlobal: int
          The global noise values for all the channels.
  _baselines: int
          Contains spike_delay number of frames of median baseline values.
  Updated by user at every frame.
  */
  if (_index_baselines < 0) {
    spikes_filtered_file.close();
    cout << "Index baselines less than 0. Terminating Spike Handler" << endl;
    exit(EXIT_FAILURE);
  }
  Parameters::aGlobal = _aGlobal;
  Parameters::baselines = _baselines;
  Parameters::index_baselines = _index_baselines;
}

void addSpike(int channel, int frame, int amplitude) {
  /*Adds a spike to the spikes_to_be_processed deque. Once the frame of the
spike to be added is
  greater than the spike_peak_duration larger than the first spike in the deque,
it will process
  all the current spikes in the deque and then attempt to add the spike again.
It will keep repeating
  this process until the spike to be added frame is smaller than the
spike_peak_duration plus the frames
of the first spike or the deque is empty.

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

  if (channel < Parameters::num_channels && channel >= 0) {
    Spike spike_to_be_added;
    spike_to_be_added.channel = channel;
    spike_to_be_added.frame = frame;
    spike_to_be_added.amplitude = amplitude;

    if (Parameters::debug) {
      cout << "storing COM cutouts " << endl;
    }
    spike_to_be_added = storeWaveformCutout(cutout_size, spike_to_be_added);
    if (Parameters::debug) {
      cout << "... done storing COM cutouts " << endl;
    }
    if (Parameters::to_localize) {
      if (Parameters::debug) {
        cout << "Storing counts..." << endl;
      }
        spike_to_be_added = storeCOMWaveformsCounts(spike_to_be_added);
        if (Parameters::debug) {
          cout << "... done storing counts!"  << endl;
        }
    }

    bool isAdded = false;
    while (!isAdded) {
      if (Parameters::spikes_to_be_processed.empty()) {
        Parameters::spikes_to_be_processed.push_back(spike_to_be_added);
        isAdded = true;
      } else {
        Spike first_spike = Parameters::spikes_to_be_processed.front();
        int first_frame = first_spike.frame;
        int current_frame = spike_to_be_added.frame;
        if (current_frame > first_frame + (Parameters::spike_peak_duration +
                                           Parameters::noise_duration)) {
          if (Parameters::to_localize) {
            try {
              if (Parameters::debug) {
                cout << "spike frame: " << spike_to_be_added.frame << endl;
              }
                ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file,
                                                    filteredsp);
            } catch (...) {
              spikes_filtered_file.close();
              cout << "Baseline matrix or its parameters entered incorrectly. "
                      "Terminating SpikeHandler."
                   << endl;
              exit(EXIT_FAILURE);
            }
          } else {
            ProcessSpikes::filterSpikes(spikes_filtered_file, filteredsp);
          }
        } else {
          Parameters::spikes_to_be_processed.push_back(spike_to_be_added);
          isAdded = true;
        }
      }
    }
  }
}
void terminateSpikeHandler() {
  // Filter any remaining spikes leftover at the end and close the spike file.
  while (Parameters::spikes_to_be_processed.size() != 0) {
    if (Parameters::to_localize) {
      ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file, filteredsp);
    } else {
      ProcessSpikes::filterSpikes(spikes_filtered_file, filteredsp);
    }
  }
  spikes_filtered_file.close();
  if (Parameters::verbose) {
    filteredsp.close();
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

void fillNeighborLayerMatrices() {
  if (Parameters::debug) {
    cout << "Filling Neighbor Layer Matrix" << endl;
  }
  int curr_channel;
  int curr_neighbor;
  float curr_dist;
  vector<tuple<int, float>> distances_neighbors;
  vector<int> inner_neighbors;
  for (int i = 0; i < Parameters::num_channels; i++) {
    curr_channel = i;
    for (int j = 0; j < Parameters::max_neighbors; j++) {
      curr_neighbor = Parameters::neighbor_matrix[curr_channel][j];
      if (curr_channel != curr_neighbor && curr_neighbor != -1) {
        curr_dist = channelsDist(curr_neighbor, curr_channel);
        distances_neighbors.push_back(make_tuple(curr_neighbor, curr_dist));
      }
    }
    if (distances_neighbors.size() != 0) {
      sort(begin(distances_neighbors), end(distances_neighbors),
           CustomLessThan());
      inner_neighbors = getInnerNeighborsRadius(distances_neighbors, i);
    }

    vector<int>::iterator it;
    it = inner_neighbors.begin();
    int curr_inner_neighbor;
    int k = 0;
    // Fill Inner neighbors matrix
    while (it != inner_neighbors.end()) {
      curr_inner_neighbor = *it;
      Parameters::inner_neighbor_matrix[i][k] = curr_inner_neighbor;
      ++k;
      ++it;
    }
    while (k < Parameters::max_neighbors) {
      Parameters::inner_neighbor_matrix[i][k] = -1;
      ++k;
    }
    // Fill outer neighbor matrix
    k = 0;
    for (int l = 0; l < Parameters::max_neighbors; l++) {
      curr_neighbor = Parameters::neighbor_matrix[i][l];
      if (curr_neighbor != -1 && curr_neighbor != i) {
        bool is_outer_neighbor = true;
        for (size_t m = 0; m < inner_neighbors.size(); m++) {
          if (Parameters::inner_neighbor_matrix[i][m] == curr_neighbor) {
            is_outer_neighbor = false;
            break;
          }
        }
        if (is_outer_neighbor) {
          Parameters::outer_neighbor_matrix[i][k] = curr_neighbor;
          ++k;
        }
      }
    }
    while (k < Parameters::max_neighbors) {
      Parameters::outer_neighbor_matrix[i][k] = -1;
      ++k;
    }
    inner_neighbors.clear();
    distances_neighbors.clear();
  }
}

vector<int>
getInnerNeighborsRadius(vector<tuple<int, float>> distances_neighbors,
                        int central_channel) {
  int curr_neighbor;
  float curr_dist;
  vector<int> inner_channels;
  vector<tuple<int, float>>::iterator it;
  inner_channels.push_back(central_channel);
  it = distances_neighbors.begin();
  while (it != distances_neighbors.end()) {
    curr_neighbor = get<0>(*it);
    curr_dist = get<1>(*it);
    if (curr_dist <= Parameters::inner_radius) {
      inner_channels.push_back(curr_neighbor);
      ++it;
    } else {
      break;
    }
  }
  return inner_channels;
}

int **createInnerNeighborMatrix() {
  int **inner_neighbor_matrix;

  inner_neighbor_matrix = new int *[Parameters::num_channels];
  for (int i = 0; i < Parameters::num_channels; i++) {
    inner_neighbor_matrix[i] = new int[Parameters::max_neighbors];
  }

  return inner_neighbor_matrix;
}

int **createOuterNeighborMatrix() {
  int **outer_neighbor_matrix;

  outer_neighbor_matrix = new int *[Parameters::num_channels];
  for (int i = 0; i < Parameters::num_channels; i++) {
    outer_neighbor_matrix[i] = new int[Parameters::max_neighbors];
  }
  return outer_neighbor_matrix;
}

Spike storeWaveformCutout(int cutout_size, Spike curr_spike) {
  /*Stores the largest waveform in the written_cutout deque of the spike

  Parameters
  ----------
  channel: int
          The channel at which the spike is detected.
  cutout_size: int
          The length of the cutout
  curr_spike: Spike
          The current detected spike

  Returns
  ----------
  curr_spike: Spike
          The current detected spike (now with the waveform stored)
  */
  int channel = curr_spike.channel;
  int32_t curr_written_reading;
  int frames_processed = Parameters::max_frames_processed * Parameters::iterations;
  for (int i = 0; i < cutout_size; i++) {
    try {
      int curr_reading_index =
          (curr_spike.frame - Parameters::cutout_start - frames_processed +
           Parameters::index_data + i) * Parameters::num_channels + channel;
      if (curr_reading_index < 0 ||
          curr_reading_index > Parameters::end_raw_data) {
        curr_written_reading = (int32_t)0;
      } else {
        curr_written_reading =
            (int32_t)Parameters::raw_data[curr_reading_index];
      }
    } catch (...) {
      spikes_filtered_file.close();
      cout << "Raw Data and it parameters entered incorrectly, could not "
              "access data. Terminating SpikeHandler."
           << endl;
      exit(EXIT_FAILURE);
    }
    curr_spike.written_cutout.push_back(curr_written_reading);
  }
  return curr_spike;
}

Spike storeCOMWaveformsCounts(Spike curr_spike) {
  /*Stores the inner neighbor waveforms in the amp_cutouts deque of the spike

  Parameters
  ----------
  channel: int
          The channel at which the spike is detected.
  curr_spike: Spike
          The current detected spike

  Returns
  ----------
  curr_spike: Spike
          The current detected spike (now with the nearest waveforms stored)
  */

  vector<int> com_cutouts;
  int *nearest_neighbor_counts;
  nearest_neighbor_counts = new int[Parameters::num_com_centers];
  for (int i = 0; i < Parameters::num_com_centers; i++) {
    nearest_neighbor_counts[i] = 0;
  }

  //Get closest channels for COM
  int channel = curr_spike.channel;
  for (int i = 0; i < Parameters::num_com_centers; i++) {
    int curr_max_channel = Parameters::inner_neighbor_matrix[channel][i];
    if (curr_max_channel == -1) {
      cout << "num_com_centers too large. Not enough inner neighbors."
           << endl;
      exit(EXIT_FAILURE);
    }
    curr_spike.largest_channels.push_back(curr_max_channel);
    int frames_processed = Parameters::max_frames_processed * Parameters::iterations;
    for (int j = 0; j < Parameters::max_neighbors; j++) {
      int curr_neighbor_channel = Parameters::inner_neighbor_matrix[curr_max_channel][j];
      // Out of inner neighbors
      if (curr_neighbor_channel != -1) {
        // Masked neighbor
        if (Parameters::masked_channels[curr_neighbor_channel] == 1) {
          nearest_neighbor_counts[i] += 1;
          //Check if noise_duration is too large in comparison to the buffer size
          int cutout_size = Parameters::cutout_start + Parameters::cutout_end + 1;
          int amp_cutout_size, cutout_start_index;
          if(Parameters::cutout_start < Parameters::noise_duration || Parameters::cutout_end < Parameters::noise_duration) {
            amp_cutout_size = cutout_size;
            cutout_start_index = Parameters::cutout_start;
          }
          else {
            amp_cutout_size = Parameters::noise_duration*2;
            cutout_start_index = Parameters::noise_duration;
          }
          for (int k = 0; k < amp_cutout_size; k++) {
            int curr_reading =
                Parameters::raw_data[(curr_spike.frame - cutout_start_index -
                                     frames_processed +
                                     Parameters::index_data + k) *
                                     Parameters::num_channels +
                                     curr_neighbor_channel];
            int curr_amp = ((curr_reading - Parameters::aGlobal) * Parameters::ASCALE -
                             Parameters::baselines[curr_neighbor_channel]
                             [Parameters::index_baselines]);
            if (curr_amp < 0) {
              com_cutouts.push_back(0);
            } else {
              com_cutouts.push_back(curr_amp);
            }
          }
        }
      }
      // Out of neighbors to add cutout for
      else {
        break;
      }
    }

  }
  curr_spike.waveformscounts = make_tuple(com_cutouts, nearest_neighbor_counts);
  return curr_spike;
}

}
