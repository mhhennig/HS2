#include "SpikeHandler.h"
#include "SpkDonline.h"

namespace SpkDonline {
Detection::Detection() {}

void Detection::InitDetection(long nFrames, int sf, int NCh,
                              long ti, long *Indices, int agl) {
  NChannels = NCh;
  tInc = ti;
  Qd = new int[NChannels];      // noise amplitude
  Qm = new int[NChannels];      // median
  Sl = new int[NChannels];      // counter for spike length
  AHP = new bool[NChannels];    // counter for repolarizing current
  Amp = new int[NChannels];     // buffers spike amplitude
  SpkArea = new int[NChannels]; // integrates over spike
  A = new int[NChannels];       // control parameter for amplifier effects
  ChInd = new int[NChannels];
  Slice = new int[NChannels];

  Sampling = sf;
  Aglobal = new int[tInc];
  for (int i = 0; i < tInc; i++)
    Aglobal[i] = agl;
  for (int i = 0; i < NChannels; i++) {
    Qd[i] = 400;
    Qm[i] = Voffset * Ascale;
    Sl[i] = 0;
    AHP[i] = false;
    Amp[i] = 0;
    A[i] = artT; // start like after an out-of-linear-regime event
    SpkArea[i] = 0;
    ChInd[i] = Indices[i];
  }

  spikeCount = 0;

}

void Detection::SetInitialParams(float * pos_mtx,
                                int * neigh_mtx, int num_channels,
                                int spike_peak_duration,
                                string file_name, int noise_duration,
                                float noise_amp_percent, float inner_radius,
                                int *_masked_channels, int max_neighbors,
                                int num_com_centers, bool to_localize,
                                int thres, int cutout_start, int cutout_end,
                                int maa, int ahpthr, int maxsl, int minsl,
                                bool decay_filtering, bool verbose) {
  // set the detection parameters
  threshold = thres;
  MinAvgAmp = maa;
  AHPthr = ahpthr;
  MaxSl = maxsl;
  MinSl = minsl;
  float **channel_positions;
  int **neighbor_matrix;
  masked_channels = _masked_channels;
  iterations  = 0;

  channel_positions = createPositionMatrix(num_channels);
  for(int i=0; i<num_channels; i++){
    channel_positions[i][0] = pos_mtx[2*i];
    channel_positions[i][1] = pos_mtx[2*i+1];
  }
  neighbor_matrix = createNeighborMatrix(num_channels, max_neighbors);
  for(int i=0; i<num_channels; i++){
    for(int j=0; j<max_neighbors; j++) {
      neighbor_matrix[i][j] = neigh_mtx[i*max_neighbors+j];
    }
  }

  Qms = new int *[num_channels];
  for (int i = 0; i < num_channels; i++) {
    Qms[i] = new int[spike_peak_duration + maxsl + 2];
    //Qms[i] = new int[spike_peak_duration + maxsl + 1];
  }

  currQmsPosition = -1;
  write_out = verbose;

  if (write_out) {
    spikes_file.open(file_name + "_detected_spikes.asc");
  }

  SpikeHandler::setInitialParameters(
      num_channels, spike_peak_duration, file_name, noise_duration,
      noise_amp_percent, inner_radius, masked_channels, channel_positions,
      neighbor_matrix, max_neighbors, num_com_centers, to_localize,
      cutout_start, cutout_end, maxsl, decay_filtering, verbose);
}

void Detection::MedianVoltage(short *vm, int tInc, int tCut) // easier to interpret, though
                                         // it takes a little longer to
                                         // run, but I'm sure there is
                                         // a faster method in C++ for
                                         // computing the median
{ // otherwise could try to take the mean also (have to ignore channels out of
  // the linear regime then) as signals are about 15% correlated
  for (int t = tCut; t < tInc + tCut; t++) { // this function wastes most of the time
    for (int i = 0; i < NChannels; i++) { // loop across channels
        if (masked_channels[i] != 0) {
          Slice[i] = vm[i + t * NChannels];   // vm [i] [t]; // contains rubbish data in sort if masked? (Rickey)
        }
    }
    sort(Slice, Slice + NChannels);
    Aglobal[t - tCut] = Slice[NChannels / 2];
  }
}

void Detection::MeanVoltage(short *vm, int tInc,
                            int tCut) // if median takes too long...
                                      // or there are only few
                                      // channnels (?)
{
  int n;
  int Vsum;

  for (int t = tCut; t < tInc + tCut; t++) {
    n = 0; // no zero division problem unless all channels are masked
    Vsum = 0;
    for (int i = 0; i < NChannels; i++) { // loop across channels
      if (masked_channels[i] != 0) {
          Vsum += (vm[i + t * NChannels]);
          n++;
      }
    }
    Aglobal[t - tCut] = Vsum / n;
  }
}

void Detection::Iterate(short *vm, long t0, int tInc, int tCut, int tCut2, int maxFramesProcessed) {
  // MeanVoltage(vm, tInc, tCut);
  int a;    // to buffer the difference between ADC counts and Qm, and basline
  int t, i; // counters
  SpikeHandler::loadRawData(vm, tCut, iterations, maxFramesProcessed, tCut, tCut2);

  ++iterations;
  // Does this need to end at tInc + tCut? (Cole+Martino) // Yes (Rickey)
  for (t = tCut; t < tInc + tCut;
       t++) { // loop over data, will be removed for an online algorithm
              // SPIKE DETECTION
    currQmsPosition += 1;
    for (i = 0; i < NChannels; i++) { // loop across channels
                                      // CHANNEL OUT OF LINEAR REGIME) {
      if (masked_channels[i] != 0) {
        a = (vm[i + t * NChannels] - Aglobal[t - tCut]) * Ascale -
            Qm[i]; // difference between ADC counts and Qm
        // UPDATE Qm and Qd
        if (a > 0) {
          if (a > Qd[i]) {
            Qm[i] += Qd[i] / Tau_m0;
            if (a < (5 * Qd[i])) {
              Qd[i]++;
            } else if ((Qd[i] > Qdmin) & (a > (6 * Qd[i]))) {
              Qd[i]--;
            }
          } else if (Qd[i] > Qdmin) { // set a minimum level for Qd
            Qd[i]--;
          }
        } else if (a < -Qd[i]) {
          Qm[i] -= Qd[i] / Tau_m0 / 2;
        }
        Qms[i][currQmsPosition % (MaxSl + Parameters::spike_peak_duration)] = Qm[i];

        a = (vm[i + t * NChannels] - Aglobal[t - tCut]) * Ascale -
            Qm[i]; // should tCut be subtracted here?? // this is correct now
        // TREATMENT OF THRESHOLD CROSSINGS
        if (Sl[i] > 0) {                     // Sl frames after peak value
          Sl[i] = (Sl[i] + 1) % (MaxSl + 1); // increment Sl[i]
          if (Sl[i] < MinSl) { // calculate area under first and second frame
                               // after spike
            SpkArea[i] += a;
          }
          // check whether it does repolarize
          else if (a < (AHPthr * Qd[i])) {
            AHP[i] = true;
          }
          // accept spikes after MaxSl frames if...
          if ((Sl[i] == MaxSl) & (AHP[i])) {
            if ((2 * SpkArea[i]) > (MinSl * MinAvgAmp * Qd[i])) {
              // increase spike count
              spikeCount += 1;

              if (t - tCut - MaxSl + 1 > 0) {
                SpikeHandler::setLocalizationParameters(
                    Aglobal[t - tCut - MaxSl + 1], Qms,
                    (currQmsPosition + 1) % (MaxSl + Parameters::spike_peak_duration));
              } else {
                SpikeHandler::setLocalizationParameters(
                    Aglobal[t - tCut], Qms,
                    (currQmsPosition + 1) % (MaxSl + Parameters::spike_peak_duration));
              }
              if (write_out) {
                spikes_file << ChInd[i] << " " << t0 - MaxSl + t - tCut + 1
                            << " " << Amp[i] << endl;
              }

              SpikeHandler::addSpike(ChInd[i], t0 - MaxSl + t - tCut + 1,
                                     Amp[i]);
            }
            Sl[i] = 0;
          }
          // check whether current ADC count is higher
          else if (Amp[i] < a) {
            Sl[i] = 1; // reset peak value
            Amp[i] = a;
            AHP[i] = false;  // reset AHP
            SpkArea[i] += a; // not resetting this one (anyway don't need to
                             // care if the spike is wide)
          }
        }
        // check for threshold crossings
        else if (a > ((threshold * Qd[i]) / 2)) {
          Sl[i] = 1;
          Amp[i] = a;
          AHP[i] = false;
          SpkArea[i] = a;
        }
      }
    }
  }

} // Iterate

void Detection::FinishDetection() // write spikes in interval after last
                                  // recalibration; close file
{
  SpikeHandler::terminateSpikeHandler();
  if (write_out) {
    spikes_file.close();
  }
  else {
    spikes_file
         << "Turn on verbose in DetectFromRaw method to get all detected spikes"
         << endl;
    spikes_file.close();
  }
}

void buildPositionsMatrix(float **_channel_positions, string positions_file_path,
                          int rows, int cols) {
  /**
Reads from a string file and fills an array that contains coordinate positions
of each channel in the probe.

Parameters
  ----------
  channel_positions: 2D int array
          A 2D array where each index has the X and Y position of every channel.
  */
  int line_count;
  float coordinate;
  string line;
  ifstream positions_file(positions_file_path);
  line_count = 0;
  string::size_type sz;

  if (positions_file.is_open()) {
    while (getline(positions_file, line)) {
      stringstream ss(line);
      string token;
      int index = 0;
      while (getline(ss, token, ',')) {
        coordinate = stof(token, &sz);
        _channel_positions[line_count][index] = coordinate;
        ++index;
      }
      ss.str(string());
      index = 0;
      ++line_count;
    }
    positions_file.close();
  }
}

void buildNeighborMatrix(int **_neighbor_matrix, string neighbors_file_path,
                         int rows, int cols) {
  /**
Reads from a string file and fills an array which contains the neighbors
of each channel. Neighbors are channels that also receive waves from the
same neural spike.

Parameters
  ----------
  neighbor_matrix: 2D int array
          A 2D array with each index representing channel numbers that
          correspond to integer array values that contain the channel
          numbers that the index channel is neighboring.
  */
  int line_count;
  int neighbor;
  string line;
  ifstream neighbor_matrix_file(neighbors_file_path);
  line_count = 0;
  string::size_type sz;

  if (neighbor_matrix_file.is_open()) {
    while (getline(neighbor_matrix_file, line)) {
      stringstream ss(line);
      string token;
      int index = 0;
      while (getline(ss, token, ',')) {
        neighbor = stoi(token, &sz);
        _neighbor_matrix[line_count][index] = neighbor;
        ++index;
      }
      while (index < cols) {
        _neighbor_matrix[line_count][index] = -1;
        ++index;
      }
      ss.str(string());
      index = 0;
      ++line_count;
    }
    neighbor_matrix_file.close();
  }
}

float **createPositionMatrix(int position_rows) {
  float **_channel_positions;

  _channel_positions = new float *[position_rows];
  for (int i = 0; i < position_rows; i++) {
    _channel_positions[i] = new float[2];
  }

  return _channel_positions;
}

int **createNeighborMatrix(int channel_rows, int channel_cols) {
  int **_neighbor_matrix;

  _neighbor_matrix = new int *[channel_rows];
  for (int i = 0; i < channel_rows; i++) {
    _neighbor_matrix[i] = new int[channel_cols];
  }

  return _neighbor_matrix;
}

int **createBaselinesMatrix(int channel_rows, int channel_cols) {
  int **_Qms;

  _Qms = new int *[channel_rows];
  for (int i = 0; i < channel_rows; i++) {
    _Qms[i] = new int[channel_cols];
  }

  return _Qms;
}
}
