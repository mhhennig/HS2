import numpy as np
import matplotlib.pyplot as plt
import math

# return list of cutouts to be ignored
def filteredSpikes(file):
	g = open('filtered_spikes', 'w')
	h = open('filtered_shapes', 'w')

	with open(file) as f:
		line_count = 0
		print_count = 0
		curr_max_line = 0
		for line in f:
			line = line.split()
			ch = int(line[0])
			frame = int(line[1])
			amp = int(line[2])
			cutout = line[3:]
			if line_count == 0:
				curr_max_ch = ch
				curr_max_frame = frame
				curr_max_amp = amp
				curr_max_cutout = ' '.join(cutout)
			else:
				if frame - curr_max_frame <= 5:
					if abs(ch - curr_max_ch) <= 7:
						if amp > curr_max_amp:
							curr_max_ch = ch
							curr_max_frame = frame
							curr_max_amp = amp
							curr_max_cutout = ' '.join(cutout)
							curr_max_line = line_count
					else:
						g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
						h.write(curr_max_cutout + '\n')
						print_count += 1
						curr_max_ch = ch
						curr_max_frame = frame
						curr_max_amp = amp
						curr_max_cutout = ' '.join(cutout)
						curr_max_line = line_count
				else:
					g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
					h.write(curr_max_cutout + '\n')
					print_count += 1
					curr_max_ch = ch
					curr_max_frame = frame
					curr_max_amp = amp
					curr_max_cutout = ' '.join(cutout)
					curr_max_line = line_count

			line_count += 1

		# include append for last line 
		g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
		h.write(curr_max_cutout + '\n')
		print_count += 1

	print "Number of filtered spikes =", print_count



def amps(cutout):
	baseline = np.median(np.asarray(cutout[:4]))
	peak = np.min(np.asarray(cutout[4:19]))
	return baseline - peak


def localisation(file, chpos, clen=26, medians=False):
	chpos = np.load(chpos)

	# coord arrays
	X = []
	Y = []

	g = open('local_spikes', 'w')
	h = open('pca_cutouts', 'w')

	with open(file) as f:
		line_count = 0
		for line in f:
			line_count += 1
			line = map(float, line.split())
			no_ch = int(len(line))/clen
			if len(line) % clen != 0:
				print "Length of line =", len(line)
				print "Number of channels =", no_ch
				print "Line number =", line_count
	 		# Create dictionary to store channel info
			# spike_dict = {}
			amps_ = []
			chpos_ = []
			max_amp = 0
			max_amp_cutout = ''
			for i in range(no_ch):
				chID = int(line[i*clen])
				cutout = np.array(line[i*clen+1:(i+1)*clen])
				cutout_amp = amps(cutout)
				if cutout_amp > max_amp:
					max_amp = cutout_amp
					max_amp_cutout = ' '.join(map(str,cutout))
				# spike_dict[chID] = [chpos[chID], cutout_amp]
				amps_.append(cutout_amp)
				chpos_.append(chpos[chID])
			if line_count == 0:
				print amps_
				print chpos_
			# Write max cutout to file
			h.write(max_amp_cutout + '\n')
			# Calculate spike medians
			amps_ = np.asarray(amps_)
			chpos_ = np.asarray(chpos_)
			if medians == True:
				# Calculate spike medians
				median = np.median(amps_)
				amps_ = np.clip(amps_-median,0.0,1.0e12)
			# Calculate centre of mass
			if np.sum(amps_) != 0:
				Bj = np.dot(amps_, chpos_)/np.sum(amps_)
				g.write(str(Bj[1]) + ' ' + str(Bj[0]) + '\n')
			# else:
				print "Broken line =", line_count
				print line


