#include "SpikeHandler.h"

int Parameters::num_channels;
int Parameters::spike_delay;
int Parameters::spike_peak_duration;
int Parameters::noise_duration;
float Parameters::noise_amp_percent;
int Parameters::max_neighbors;
int** Parameters::neighbor_matrix;
int** Parameters::channel_positions;
int** Parameters::inner_neighbor_matrix;
int** Parameters::outer_neighbor_matrix;
int Parameters::aGlobal;
bool Parameters::to_localize;
bool Parameters::verbose;
int** Parameters::baselines;
int Parameters::cutout_start;
int Parameters::cutout_end;
int Parameters::index_data;
int Parameters::index_baselines;
int Parameters::iterations;
int Parameters::frames;
int Parameters::maxsl;
int Parameters::end_raw_data;
short* Parameters::raw_data;
int* Parameters::masked_channels;
int Parameters::event_number;
bool Parameters::debug;
float Parameters::inner_radius;

deque<Spike> Parameters::spikes_to_be_processed;
std::ofstream filteredsp;
std::ofstream spikes_filtered_file;


void setInitialParameters(int _num_channels, int _spike_delay, int _spike_peak_duration, string file_name, \
						  int _noise_duration, float _noise_amp_percent, float _inner_radius, int* _masked_channels, int** _channel_positions, int** _neighbor_matrix, \
						  int _max_neighbors, bool _to_localize = false, int _cutout_start= 10, int _cutout_end=20, int _maxsl = 0, bool _verbose = false)
{
	/*This sets all the initial parameters needed to run the filtering algorithm.

	Parameters
	----------
	_num_channels: int
		Number of channels on the probe
	_spike_delay: int
		The number of frames back a spike occurred after it was detected (where the beginning of the spike was).
	file_name: string
        The name of the file that the processed spikes will be written to in binary.
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
    _masked_channels: 1D int array
    	The channels that are masked for the detection (AKA we get input from them)
	_channel_positions: 2D int array
		Indexed by the channel number starting at 0 and going up to num_channels - 1. Each
		index contains pointer to another array which contains X and Y position of the channel. User creates
		this before calling SpikeHandler.
	_neighbor_matrix: 2D int array
		Indexed by the channel number starting at 0 and going up to num_channels - 1. Each
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
	if(_num_channels < 0) {
		cout << "Number of channels given incorrectly. Terminating Spike Handler" << endl;
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
    if(_inner_radius < 0) {
		cout << "Inner Radius less than 0. Terminating Spike Handler" << endl;
		exit(EXIT_FAILURE);
	}

	Parameters::num_channels = _num_channels;
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
    Parameters::masked_channels = _masked_channels;
    Parameters::inner_radius = _inner_radius;
    Parameters::event_number = 0;
    Parameters::debug = false;
    Parameters::verbose = _verbose;

    Parameters::inner_neighbor_matrix = createInnerNeighborMatrix();
    Parameters::outer_neighbor_matrix = createOuterNeighborMatrix();
    fillNeighborLayerMatrices();
    if(Parameters::verbose) {
        for(int i=0; i<Parameters::num_channels; i++) {
            cout << "Channel: " << i << endl;
            cout << "Inner Neighbors: ";
            for(int j=0; j < Parameters::max_neighbors - 1; j++) {
         	      cout << Parameters::inner_neighbor_matrix[i][j]  << "  ";
            }
            cout << endl;
            cout << "Outer Neighbors: ";
            for(int k=0; k < Parameters::max_neighbors - 1; k++) {
        	       cout <<  Parameters::outer_neighbor_matrix[i][k]  << "  ";
            }
    	    cout << endl;
        }
    }

	spikes_filtered_file.open(file_name, ios::binary);
    filteredsp.open("Filtered Spikes");

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
	this process until the spike to be added frame is smaller than the spike_peak_duration plus the frames
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
	int amp_cutout_size = Parameters::spike_delay*2 + 1;
	int frames_processed = Parameters::frames*Parameters::iterations;

	if(channel < Parameters::num_channels && channel >= 0) {
		int curr_neighbor_channel, curr_reading;
		int32_t curr_written_reading;
		Spike spike_to_be_added;
		spike_to_be_added.channel = channel;
		spike_to_be_added.frame = frame;
		spike_to_be_added.amplitude = amplitude;
		int ASCALE = -64;

		for(int i = 0; i < cutout_size; i++) {
			try {
				int curr_reading_index = (frame - Parameters::cutout_start - frames_processed + Parameters::index_data + i)*Parameters::num_channels + channel;
				if(curr_reading_index < 0 || curr_reading_index > Parameters::end_raw_data) {
					curr_written_reading = (int32_t) 0;
				}
				else {
					curr_written_reading = (int32_t) Parameters::raw_data[curr_reading_index];
				}
			} catch (...) {
				spikes_filtered_file.close();
				cout << "Raw Data and it parameters entered incorrectly, could not access data. Terminating SpikeHandler." << endl;
				exit(EXIT_FAILURE);
			}
			spike_to_be_added.written_cutout.push_back(curr_written_reading);
        }

		for(int i = 0; i < Parameters::max_neighbors - 1; i++) {
			try {
				curr_neighbor_channel = Parameters::inner_neighbor_matrix[channel][i];
			} catch (...) {
				spikes_filtered_file.close();
				cout << "Neighbor matrix improperly created. Terminating SpikeHandler" << endl;
				exit(EXIT_FAILURE);
			}
            //Out of inner neighbors
			if(curr_neighbor_channel != -1) {
                //Masked neighbor
                if(Parameters::masked_channels[curr_neighbor_channel] == 1) {
    				for(int j = 0; j < amp_cutout_size; j++) {
    					try {
      						curr_reading = Parameters::raw_data[(frame - Parameters::spike_delay - frames_processed + Parameters::index_data + j)*Parameters::num_channels + curr_neighbor_channel];
    					} catch (...) {
    						spikes_filtered_file.close();
    						cout << "Raw Data and it parameters entered incorrectly, could not access data. Terminating SpikeHandler." << endl;
    						exit(EXIT_FAILURE);
    					}

    					int curr_amp = ((curr_reading - Parameters::aGlobal) * ASCALE - Parameters::baselines[curr_neighbor_channel][Parameters::index_baselines]);
                        if(curr_amp > 100000) {
                            cout << "CURR AMP TOO BIG: " << curr_amp << endl;
                            cout << "CURR READING: " << curr_reading << endl;
                            cout << "Baseline: " << Parameters::baselines[curr_neighbor_channel][Parameters::index_baselines] << endl;
                            cout << "Global: " << Parameters::aGlobal << endl;
                            cout << "Scale: " << ASCALE << endl;
                        }
                        if(curr_amp < 0) {
                            spike_to_be_added.amp_cutouts.push_back(0);
                        } else {
                            spike_to_be_added.amp_cutouts.push_back(curr_amp);
                        }
    				}
                }
			}
			//Out of neighbors to add cutout for
			else {
				break;
			}
		}
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
                            if(Parameters::debug && spike_to_be_added.frame > 120) {
                                ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file, filteredsp);
                                spikes_filtered_file.close();
                                filteredsp.close();
    							cout << "Baseline matrix or its parameters entered incorrectly. Terminating SpikeHandler." << endl;
    							exit(EXIT_FAILURE);
                            }
                            else {
                                ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file, filteredsp);;
                            }
						} catch (...) {
							spikes_filtered_file.close();
							cout << "Baseline matrix or its parameters entered incorrectly. Terminating SpikeHandler." << endl;
							exit(EXIT_FAILURE);
						}
					}
					else {
						ProcessSpikes::filterSpikes(spikes_filtered_file,filteredsp);
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
			ProcessSpikes::filterLocalizeSpikes(spikes_filtered_file, filteredsp);
		}
		else {
			ProcessSpikes::filterSpikes(spikes_filtered_file, filteredsp);
		}
	}
	spikes_filtered_file.close();
    if(!Parameters::verbose) {
        filteredsp << "Turn on verbose in DetectFromRaw method to get all filtered spikes" << endl;
    }
    filteredsp.close();
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
    if(Parameters::debug) {
        cout << "Filling Neighbor Layer Matrix" << endl;
    }
    int curr_channel;
    int curr_neighbor;
    float curr_dist;
    vector<tuple<int, float>> distances_neighbors;
    vector<int> inner_neighbors;
    for(int i = 0; i < Parameters::num_channels; i++) {
        curr_channel = i;
        for(int j = 0; j < Parameters::max_neighbors; j++) {
            curr_neighbor = Parameters::neighbor_matrix[curr_channel][j];
            if(curr_channel != curr_neighbor && curr_neighbor != -1) {
                curr_dist = channelsDist(curr_neighbor, curr_channel);
                distances_neighbors.push_back(make_tuple(curr_neighbor, curr_dist));
            }
        }
        if(distances_neighbors.size() != 0) {
            sort(begin(distances_neighbors), end(distances_neighbors), CustomLessThan());
            if(Parameters::debug) {
                cout << "Found Distances" << endl;
            }
            //inner_neighbors = getInnerNeighborsBounding(distances_neighbors, i);
            inner_neighbors = getInnerNeighborsRadius(distances_neighbors, i);
        }

        if(Parameters::debug) {
            cout << "Got inner neighbors" << endl;
            cout << i << endl;
        }

        vector<int>::iterator it;
        it = inner_neighbors.begin();
        int curr_inner_neighbor;
        int k = 0;
        //Fill Inner neighbors matrix
        while(it != inner_neighbors.end())
        {
            curr_inner_neighbor = *it;
            Parameters::inner_neighbor_matrix[i][k] = curr_inner_neighbor;
            ++k;
            ++it;
        }
        while(k < Parameters::max_neighbors - 1) {
            Parameters::inner_neighbor_matrix[i][k] = -1;
            ++k;
        }
        if(Parameters::debug) {
            cout << "Filling Inner Neighbors" << endl;
            cout << i << endl;
        }
        //Fill outer neighbor matrix
        k = 0;
        for(int l = 0; l < Parameters::max_neighbors; l++) {
            curr_neighbor = Parameters::neighbor_matrix[i][l];
            if(curr_neighbor != -1 && curr_neighbor != i) {
                bool is_outer_neighbor = true;
                for(size_t m = 0; m < inner_neighbors.size(); m++) {
                    if(Parameters::inner_neighbor_matrix[i][m] == curr_neighbor) {
                        is_outer_neighbor = false;
                        break;
                    }
                }
                if(is_outer_neighbor) {
                    Parameters::outer_neighbor_matrix[i][k] = curr_neighbor;
                    ++k;
                }
            }
        }
        while(k < Parameters::max_neighbors - 1) {
            Parameters::outer_neighbor_matrix[i][k] = -1;
            ++k;
        }
        inner_neighbors.clear();
        distances_neighbors.clear();
    }
}

vector<int> getInnerNeighborsRadius(vector<tuple<int, float>> distances_neighbors, int central_channel) {
    int curr_neighbor;
    float curr_dist;
    vector<int> inner_channels;
    vector<tuple<int, float>>::iterator it;
    it = distances_neighbors.begin();
    while(it != distances_neighbors.end())
    {
        curr_neighbor = get<0>(*it);
        curr_dist = get<1>(*it);
        if(curr_dist <= Parameters::inner_radius) {
            inner_channels.push_back(curr_neighbor);
            ++it;
        }
        else {
            break;
        }

    }
    return inner_channels;
}

vector<int> getInnerNeighborsBounding(vector<tuple<int, float>> distances_neighbors, int central_channel) {
    vector<int> inner_neighbors;
    vector<Point> boundary_points;
    vector<Line> boundary_lines;
    vector<tuple<int, float>>::iterator it;
    it = distances_neighbors.begin();
    int curr_neighbor;
    int central_x = Parameters::channel_positions[central_channel][0];
    int central_y = Parameters::channel_positions[central_channel][1];
    Point central_point = {central_x, central_y, central_channel};
    Point curr_point;
    Point prev_point;
    int curr_x;
    int curr_y;
    while(it != distances_neighbors.end())
    {
        curr_neighbor = get<0>(*it);
        curr_x = Parameters::channel_positions[curr_neighbor][0];
        curr_y = Parameters::channel_positions[curr_neighbor][1];
        curr_point = {curr_x, curr_y, curr_neighbor};


        if(boundary_points.size() == 0) {
            boundary_points.push_back(curr_point);
            ++it;
        }
        else if(boundary_points.size() == 1) {
            prev_point = boundary_points.front();
            Line newLine = createLine(curr_point, prev_point);
            boundary_lines.push_back(newLine);
            boundary_points.push_back(curr_point);
            ++it;
        }
        else {
            if(acceptAsBoundaryPoint(curr_point, central_point, boundary_lines)) {
                boundary_points.push_back(curr_point);
                createBoundaryLines(curr_point, boundary_points, boundary_lines);
                ++it;
            } else {
                ++it;
            }
        }
    }
    vector<int> inner_channels;
    vector<Point>::iterator it2;
    it2 = boundary_points.begin();
    Point curr_bp;
    int curr_channel;
    while(it2 != boundary_points.end())
    {
        curr_bp = *it2;
        curr_channel = curr_bp.channel;
        inner_channels.push_back(curr_channel);
        ++it2;
    }
    return inner_channels;
}

void createBoundaryLines(Point curr_point, vector<Point> &boundary_points, vector<Line> &boundary_lines) {
    Point prev_bp;
    for(size_t i = 0; i < boundary_points.size(); i++) {
        prev_bp = boundary_points.at(i);
        Line newLine = createLine(curr_point, prev_bp);
        boundary_lines.push_back(newLine);
    }
}

bool acceptAsBoundaryPoint(Point curr_point, Point central_point, vector<Line> &boundary_lines) {
    bool accept_bp = true;
    Line curr_bl;
    float line_dist_from_curr_point;
    float line_dist_from_central_point;
    vector<Line>::iterator it;
    it = boundary_lines.begin();
    while(it != boundary_lines.end())
    {
        curr_bl = *it;
        line_dist_from_curr_point = curr_bl.a*curr_point.x + curr_bl.b*curr_point.y + curr_bl.c;
        line_dist_from_central_point = curr_bl.a*central_point.x + curr_bl.b*central_point.y + curr_bl.c;
        //On same side of boundary line
        if(line_dist_from_curr_point < 0 && line_dist_from_central_point < 0) {
            accept_bp = true;
        }
        //On same side of boundary line
        else if(line_dist_from_curr_point > 0 && line_dist_from_central_point > 0) {
            accept_bp = true;
        }
        else if(line_dist_from_central_point == 0) {
            accept_bp = true;
        }
        else if(line_dist_from_curr_point == 0) {
            accept_bp = true;
        }
        //On opposite sides of boundary line
        else {
            accept_bp = false;
            break;
        }
        ++it;
    }
    return accept_bp;
}

float distBetweenPoints(Point p1, Point p2) {
    float x_displacement;
    float y_displacement;
    float dist;
    x_displacement = p1.x - p2.x;
    y_displacement = p1.y - p2.y;
    dist = sqrt(pow(x_displacement, 2) + pow(y_displacement, 2));

    return dist;
}

Line createLine(Point p1, Point p2) {
    Line newLine;
    float dx, dy;
    float slope;
    float y_intercept;
    dx = p1.x - p2.x;
    dy = p1.y - p2.y;
    if(dx == 0) {
        newLine.p1 = p1;
        newLine.p2 = p2;
        newLine.a = 1;
        newLine.b = 0;
        newLine.c = -p1.x;
    }
    else if(dy == 0) {
        newLine.p1 = p1;
        newLine.p2 = p2;
        newLine.a = 0;
        newLine.b = 1;
        newLine.c = -p1.y;
    }
    else {
        slope = dy/dx;
        y_intercept = p1.y - slope*p1.x;
        newLine.p1 = p1;
        newLine.p2 = p2;
        newLine.a = -slope;
        newLine.b = 1;
        newLine.c = -y_intercept;
    }
    return newLine;
}



int** createInnerNeighborMatrix() {
    int ** inner_neighbor_matrix;

    inner_neighbor_matrix = new int*[Parameters::num_channels];
    for (int i = 0; i < Parameters::num_channels; i++) {
            inner_neighbor_matrix[i] = new int[Parameters::max_neighbors - 1];
    }

    return inner_neighbor_matrix;
}

int** createOuterNeighborMatrix() {
    int ** outer_neighbor_matrix;

    outer_neighbor_matrix = new int*[Parameters::num_channels];
    for (int i = 0; i < Parameters::num_channels; i++) {
            outer_neighbor_matrix[i] = new int[Parameters::max_neighbors - 1];
    }
    return outer_neighbor_matrix;
}
