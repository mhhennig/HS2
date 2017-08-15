import numpy as np
import matplotlib.pyplot as plt
import math
from datetime import datetime
from bitarray import bitarray

# return list of cutouts to be ignored
def filteredSpikes(spikes, shapes, count):
	startTime = datetime.now()

	g = open('filtered_spikes', 'w')
	h = open('filtered_shapes', 'w')

	with open(count) as f:
		for line in f:
			spikeCount = int(line)

	ba = bitarray(spikeCount)
	ba.setall(False) 

	with open(spikes) as f:
		line_count = 0
		print_count = 0
		curr_max_line = 0
		for line in f:
			line = line.split()
			ch = int(line[0])
			frame = int(line[1])
			amp = int(line[2])
			if line_count == 0:
				curr_max_ch = ch
				curr_max_frame = frame
				curr_max_amp = amp
			else:
				if frame - curr_max_frame <= 5:
					if abs(ch - curr_max_ch) <= 7:
						if amp > curr_max_amp:
							curr_max_ch = ch
							curr_max_frame = frame
							curr_max_amp = amp
							curr_max_line = line_count
					else:
						g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
						ba[curr_max_line] = 1
						print_count += 1
						curr_max_ch = ch
						curr_max_frame = frame
						curr_max_amp = amp
						curr_max_line = line_count
				else:
					g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
					ba[curr_max_line] = 1
					print_count += 1
					curr_max_ch = ch
					curr_max_frame = frame
					curr_max_amp = amp
					curr_max_line = line_count

			line_count += 1

		# include append for last line
		g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
		ba[curr_max_line] = 1
		print_count += 1

	#loop through spikes
	with open(shapes) as f:
		line_count = 0
		for line in f:
			if ba[line_count] == 1:
				h.write(line)
			line_count += 1

	endTime = datetime.now()

	print("Number of filtered spikes = " + str(print_count))
	print("Time taken = " + str(endTime - startTime))



def amplitude(cutout, variability):
	z = np.hstack((cutout[:15], cutout[-30:]))
	a = (np.median(z) - np.min(cutout[25:35]))/variability
	if a > 0:
		return a
	else:
		return 0


def localisation(file, chpos, medians=True):
	startTime = datetime.now()

	chpos = np.load(chpos)

	clen = 99

	# coord arrays
	X = []
	Y = []

	g = open('local_spikes', 'w')
	h = open('pca_cutouts', 'w')

	with open(file) as f:
		line_count = 0
		count2 = 0
		count_emptystring = 0
		reference_count = 0
		zero_amps = 0
		for line in f:
			line_count += 1
			line = map(float, line.split())
			baseline = line[0]
			cutouts = line[1:]
			no_ch = int(len(cutouts))/clen
			if no_ch == 0:
				reference_count += 1
				continue
			if len(cutouts) % clen != 0:
				print("Length of line =", len(cutouts))
				print("Number of channels =", no_ch)
				print("Line number =", line_count)
				break
	 		# Create arrays to store channel info
			amps_ = []
			chpos_ = []
			max_amp = 0
			max_amp_cutout = ''
			for i in range(no_ch):
				chID = int(cutouts[i*clen])
				variability = float(cutouts[i*clen+1])
				cutout = np.array(cutouts[i*clen+2:(i+1)*clen])
				cutout_amp = amplitude(cutout, variability)
				if cutout_amp >= max_amp:
					max_amp = cutout_amp
					# max_amp_cutout = ' '.join(map(str,cutout[20:51]-baseline))
					max_amp_cutout = ' '.join(map(str,cutout[20:80]-baseline))
				amps_.append(cutout_amp)
				chpos_.append(chpos[chID])
			# Write max cutout to file
			if max_amp_cutout == '':
				# print "max_amp_cutout == '':" + str(amps_)
				count_emptystring += 1
				continue
			# Calculate spike medians
			amps_ = np.asarray(amps_)
			chpos_ = np.asarray(chpos_)
			if medians == True:
				# Calculate spike medians
				median = np.median(amps_)
				amps_1 = np.clip(amps_-median,0.0,1.0e12)
			# Count number of points that fall on grid (mean position of 2 channels)
			if np.count_nonzero(amps_1) <= 2:
				count2 += 1
			# Calculate centre of mass
			if np.sum(amps_1) != 0:
				Bj = np.dot(amps_1, chpos_)/np.sum(amps_1)
				g.write(str(Bj[1]) + ' ' + str(Bj[0]) + '\n')
				h.write(max_amp_cutout + '\n')
			else:
				print("Broken line =", line_count)
				print(amps_)
				zero_amps += 1

		endTime = datetime.now()

		print("count2 = " + str(count2))
		print("empty string = " + str(count_emptystring))
		print("reference count = " + str(reference_count))
		print("zero amps =" + str(zero_amps))
		print("Time taken = " + str(endTime - startTime))
