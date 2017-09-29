#include "SpkDonline.h"
#include "SpikeHandler.h"

namespace SpkDonline {
Detection::Detection() {}

void Detection::InitDetection(long nFrames, double nSec, int sf, int NCh, long ti, long *Indices, int agl, int tpref, int tpostf) {
  NChannels = NCh;
  tInc = ti;
  Qd = new int[NChannels];      // noise amplitude
  Qm = new int[NChannels];      // median
  int** Qms;    // Stores 5 frames of medians ending with current median
  Sl = new int[NChannels];      // counter for spike length
  AHP = new bool[NChannels];    // counter for repolarizing current
  Amp = new int[NChannels];     // buffers spike amplitude
  SpkArea = new int[NChannels]; // integrates over spike
  A = new int[NChannels];       // control parameter for amplifier effects
  ChInd = new int[NChannels];
  Slice = new int[NChannels];

  // Qdmean = new int[NChannels];
  // MaxSl = 8; //sf / 1000 + 1;
  // MinSl = 3; //sf / 3000 + 2;
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
    // Qdmean[i] = 0;
  }

  spikeCount = 0;
  // frameCount = 0;

  fpre = tpref;
  fpost = tpostf;



  // 10 neighbours
  // for (int i = 0; i < NChannels-1; i++) {
  //   if (i % 2 == 0) {
  //     for (int j = 0; j < 10; j++) {
  //       if (i - 4 + j >= 0 && i - 4 + j < NChannels-1) {
  //         ChInd10[i][j] = i - 4 + j;
  //       }
  //     }
  //   } else if (i % 2 == 1) {
  //     for (int j = 0; j < 10; j++) {
  //       if (i - 5 + j >= 0 && i - 5 + j < NChannels-1) {
  //         ChInd10[i][j] = i - 5 + j;
  //       }
  //     }
  //   }
  // }
  // Print out ChInd10
  // std::cout << "\n";
  // for (int i =0; i < NChannels - 1; i++) {
  //   for (int j=0; j < 10; j++) {
  //     std::cout << ChInd10[i*10 + j] << " ";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << "init";
}

void Detection::SetInitialParams(int num_channels, int num_recording_channels, int spike_delay, int spike_peak_duration, int noise_duration, \
                         		 int noise_amp, int max_neighbors, int cutout_length, bool to_localize, int thres, int maa, int ahpthr, int maxsl,
                                 int minsl) {
  // set the detection parameters
  threshold = thres;
  MinAvgAmp = maa;
  AHPthr = ahpthr;
  MaxSl = maxsl;
  MinSl = minsl;
  int** channel_positions;
  int**neighbor_matrix;
  channel_positions = createPositionMatrix(num_recording_channels);
  neighbor_matrix = createNeighborMatrix(num_recording_channels, max_neighbors);
  buildPositionsMatrix(channel_positions, "positions", num_recording_channels, 2);
  buildNeighborMatrix(neighbor_matrix, "neighbormatrix", num_recording_channels, max_neighbors);
  Qms = createBaselinesMatrix(num_channels, spike_peak_duration + 1);

  setInitialParameters(num_channels, num_recording_channels, spike_delay, spike_peak_duration, noise_duration, \
									   noise_amp, channel_positions, neighbor_matrix, max_neighbors, to_localize, cutout_length);
}

// don't know how to make multiple threads write to the same file,
// maybe you'll have to buffer values during the detection (make Iterate a list
// instead of a void)
// and then write spikes for each block
void Detection::openSpikeFile(const char *name) {
  cout << "# Writing to: " << name << "\n";
  //w.open(name);
  // fs = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
  // w = new StreamWriter(fs);
}

void Detection::openFiles(const char *spikes){
  //w.open(spikes);
}

void Detection::MedianVoltage(short *vm) // easier to interpret, though
                                                  // it takes a little longer to
                                                  // run, but I'm sure there is
                                                  // a faster method in C++ for
                                                  // computing the median
{ // otherwise could try to take the mean also (have to ignore channels out of
  // the linear regime then) as signals are about 15% correlated
  for (int t = 0; t < tInc; t++) { // this function wastes most of the time
    for (int i = 0; i < NChannels; i++) { // loop across channels
      Slice[i] = vm[i + t*NChannels];        // vm [i] [t];
    }
    sort(Slice, Slice + sizeof Slice / sizeof Slice[0]);
    Aglobal[t] = Slice[NChannels / 2];
  }
}

void Detection::MeanVoltage(short *vm, int tInc, int tCut) // if median takes too long...
                                                // or there are only few
                                                // channnels (?)
{
  int n;
  int Vsum;

  for (int t = tCut; t < tInc + tCut; t++) {
    n = 1; // constant offset doesn't matter, avoid zero division
    Vsum = 0;
    for (int i = 0; i < NChannels; i++) { // loop across channels
      // if (((vm[i * tInc + t] + 4) % 4096) > 10) {
      if (i + t*NChannels > (tInc + tCut)*NChannels) {
        cout << "line 125" << "\n";
      }
      Vsum += (vm[i + t*NChannels]);
      n++;
      // }
    }
    if (t-tCut > tInc) {
      cout << "line 133" << "\n";
    }
    Aglobal[t-tCut] = Vsum / n;
  }
}

void Detection::Iterate(short *vm, long t0, int tInc, int tCut, int tCut2) {
  // MeanVoltage(vm, tInc, tCut);
  int a, b=0; // to buffer the difference between ADC counts and Qm, and basline
  int CurrNghbr;
  // std::cout << NChannels << " " << t0 << " " << tInc << "\n";
  // std::cout.flush();
  int currQmsPosition = -1;
  loadRawData(vm, iterations, 100000, tCut);
  ++iterations;
  //cout << "iterations: " << iterations << '\n';
  for (int t = tCut; t < tInc + tCut;
       t++) { // loop over data, will be removed for an online algorithm
              // SPIKE DETECTION
    // frameCount += 1;
    // std::cout << "\n";
    // for (int x = 0; x < NChannels; x++) {
    //   std::cout << Qd[x] << " ";
    // }
    // std::cout << "\n"
    currQmsPosition += 1; 
    for (int i = 0; i < NChannels; i++) { // loop across channels
                                          // CHANNEL OUT OF LINEAR REGIME
      // if (((vm[i + t*NChannels] + 4) % NChannels) < 10) {
      //   if (A[i] <
      //       artT) { // reset only when it starts leaving the linear regime
      //     Sl[i] = 0;
      //     // std::cout << "yes it enters" << "\n"
      //     A[i] = artT;
      //   }
      // }
      // DEFAULT OPERATIONS
      // else if (A[i] == 0) {
      // if (A[i] == 0) {
        if (t-tCut >= tInc) {
          cout << "line 154: referencing index too large" << "\n";
        }
        a = (vm[i + t*NChannels] - Aglobal[t-tCut]) * Ascale - // should tCut be subtracted here??
        // a = (vm[i + t*NChannels] - Aglobal[t]) * Ascale -
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
        Qms[i][currQmsPosition % 6] = Qm[i];
        a = (vm[i + t*NChannels] - Aglobal[t-tCut]) * Ascale - Qm[i]; // should tCut be subtracted here??
        // TREATMENT OF THRESHOLD CROSSINGS
        if (Sl[i] > 0) { // Sl frames after peak value
          //std::cout << "*";
          // default
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

              // Write spikes to file
              int correctBaseline = Qms[i][(currQmsPosition - 5) % 6];
              //w << ChInd[i] << " " << t0 + t - MaxSl - tCut + 1 << " "
              //  << -Amp[i] * Ascale / Qd[i] <<   "\n";
              setLocalizationParameters(Aglobal[t-tCut], Qms);
			        addSpike(ChInd[i], t0 + t - tCut + 1, -Amp[i] * Ascale / Qd[i]);


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
            if (t-tCut >= tInc) {
            cout << "line 223: referencing index too large" << "\n";
            }
            b = Aglobal[t - tCut];// Qm[i]; // Again, should tCut be subtracted here?
            // b = Aglobal[t];
          }
        }
        // check for threshold crossings
        else if (a > ((threshold * Qd[i]) / 2)) {
          Sl[i] = 1;
          Amp[i] = a;
          AHP[i] = false;
          SpkArea[i] = a;
        }
      // }
      // // AFTER CHANNEL WAS OUT OF LINEAR REGIME
      // else {
      //   Qm[i] = (2 * Qm[i] + (vm[i + t*NChannels] - Aglobal[t]) * Ascale +
      //            2 * Qd[i]) /
      //           3; // update Qm
      //   A[i]--;
      // }
      // for (int v = 0; v <NChannels; v++) {
      //   Qdmean[v] = (frameCount * Qdmean[v] + Qd[v])/(frameCount + 1);
      // }
    }
  }
  // for (int i = 0; i < NChannels; i++) { // reset params after each chunk
  //   Qd[i] = 400;
  //   Qm[i] = Voffset * Ascale;
  //   Sl[i] = 0;
  //   AHP[i] = false;
  //   Amp[i] = 0;
  //   A[i] = artT; // start like after an out-of-linear-regime event
  //   SpkArea[i] = 0;
  // }
} // Iterate

void Detection::FinishDetection() // write spikes in interval after last
                                  // recalibration; close file
{
	terminateSpikeHandler();
 // w.close();
  //wCount.open("count");
  //wCount << spikeCount;
  //wCount.close();
  // wVar.open("variability");
  // for (int i = 0; i < NChannels; i++) {
  //   wVar << Qdmean[i] << " ";
  // }
  // wVar.close();
}

void buildPositionsMatrix(int** _channel_positions, string positions_file_path, int rows, int cols)
{
	/**
    Reads from a string file and fills an array that contains coordinate positions 
    of each channel in the probe.

    Parameters
	----------
	channel_positions: 2D int array
		A 2D array where each index has the X and Y position of every channel.
	*/
	int line_count;
	int coordinate;
	string line;
	ifstream positions_file (positions_file_path);
	line_count = 0; 
	string::size_type sz;

	if (positions_file.is_open())
	{
    	while ( getline (positions_file,line) )
    	{
    		stringstream ss(line);
    		string token;
    		int index = 0;
    		while (getline(ss,token, ','))
    		{
    			coordinate = stoi(token,&sz);
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

void buildNeighborMatrix(int** _neighbor_matrix, string neighbors_file_path, int rows, int cols)
{
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
	ifstream neighbor_matrix_file (neighbors_file_path);
	line_count = 0; 
	string::size_type sz;

	if (neighbor_matrix_file.is_open())
	{
    	while ( getline (neighbor_matrix_file,line) )
    	{
    		stringstream ss(line);
    		string token;
    		int index = 0;
    		while (getline(ss,token, ','))
    		{
    			neighbor = stoi(token,&sz);
    			_neighbor_matrix[line_count][index] = neighbor;
    			++index;
			}
			while(index < cols) {
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

int** createPositionMatrix(int position_rows) {
        int **_channel_positions;

        _channel_positions = new int*[position_rows];
        for (int i = 0; i < position_rows; i++) {
                _channel_positions[i] = new int[2];
        }

        return _channel_positions;
}

int** createNeighborMatrix(int channel_rows, int channel_cols) {
        int **_neighbor_matrix;

        _neighbor_matrix = new int*[channel_rows];
        for (int i = 0; i < channel_rows; i++) {
                _neighbor_matrix[i] = new int[channel_cols];
        }

        return _neighbor_matrix;
}

int** createBaselinesMatrix(int channel_rows, int channel_cols) {
        int **_Qms;

        _Qms = new int*[channel_rows];
        for (int i = 0; i < channel_rows; i++) {
                _Qms[i] = new int[channel_cols];
        }

        return _Qms;
}
}
