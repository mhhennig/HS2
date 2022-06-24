#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

namespace SpkDonline {
class Detection {
  int NChannels; // number of channels; is set when reading the data
  int *ChInd;    // indices for parallelization
  // Variables for variance and mean
  int *Qd; // noise amplitude
  int *Qm; // median
  int **Qms; //stores spike_delay + MaxSl baseline values;
  int* masked_channels; //stores all masked channels as 0 and regular channels as 1
  int iterations = 0;
  // Variables for the spike detection
  int *Sl;      // counter for spike length
  bool *AHP;    // counter for repolarizing current
  int *Amp;     // buffers spike amplitude
  int *SpkArea; // integrates over spike
  // Parameters for variance and mean updates
  const int Tau_m0 = 4;  // timescale for updating Qm (increment is Qd/Tau_m)
  const int Qdmin = 200; // set minimum value of Qd
  // Parameters for spike detection
  int threshold; // threshold to detect spikes >11 is likely to be real
                     // spikes, but can and should be sorted afterwards
  int AHPthr; // signal should go below that threshold within MaxSl-MinSl frames
  int MaxSl;     // dead time in frames after peak, used for further testing
  int MinAvgAmp; // minimal avg. amplitude of peak (in units of Qd)
  int MinSl;     // length considered for determining avg. spike amplitude
  // Parameters for
  long tInc; // 100, increment for reading data, has to be changed in main
             // program as well
  const int Ascale = -64; // factor to multiply to raw traces to increase
                          // resolution; definition of ADC counts had been
                          // changed!
  const int Voffset = 0;  // mean ADC counts, as initial value for Qm
  // Parameters for recalibration events and artefact handling
  const int artT = 10; // to use after artefacts; to update Qm for 10 frames
  int *A;              // control parameter for amplifier effects
  // Files to save the spikes etc.
  int Sampling;
  int *Aglobal;
  int *Slice;
  int a; // buffer for Iterate()

  int spikeCount;
  int currQmsPosition;
  bool debugging = true;
  bool write_out = false;
  std::ofstream spikes_file;


public:
  Detection();
  ~Detection();
  void InitDetection(long nFrames, int sf, int NCh, long ti, long int *Indices, int agl);
  void SetInitialParams(float * pos_mtx, int * neigh_mtx, int num_channels, int spike_peak_duration,
                        string file_name, int noise_duration, float noise_amp_percent, float inner_radius, int* _masked_channels, int max_neighbors,
                        int num_com_centers, bool to_localize, int thres, int cutout_start, int cutout_end, int maa, int ahpthr, int maxsl, int minsl,
                        bool decay_filtering, bool verbose);
  void MedianVoltage(short *vm, int tInc, int tCut);
  void MeanVoltage(short *vm, int tInc, int tCut);
  void Iterate(short *vm, long t0, int tInc, int tCut, int tCut2, int maxFramesProcessed);
  void FinishDetection();
};
  void buildPositionsMatrix(float** _channel_positions, string positions_file_path, int rows, int cols);
  void buildNeighborMatrix(int** _neighbor_matrix, string neighbors_file_path, int rows, int cols);
  float** createPositionMatrix(int position_rows);
  int** createNeighborMatrix(int channel_rows, int channel_cols);
int** createBaselinesMatrix(int channel_rows, int channel_cols);
};
