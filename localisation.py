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
			baseline = int(line[3])
			cutout = list(map(int,line[4:]))
			if line_count == 0:
				curr_max_ch = ch
				curr_max_frame = frame
				curr_max_amp = amp
				curr_max_cutout = ' '.join(map(str,cutout))
			else:
				if frame - curr_max_frame <= 5:
					if abs(ch - curr_max_ch) <= 7:
						if amp > curr_max_amp:
							curr_max_ch = ch
							curr_max_frame = frame
							curr_max_amp = amp
							curr_max_cutout = ' '.join(map(str,cutout))
							curr_max_line = line_count
					else:
						g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
						h.write(str(baseline) + ' ' + curr_max_cutout + '\n')
						print_count += 1
						curr_max_ch = ch
						curr_max_frame = frame
						curr_max_amp = amp
						curr_max_cutout = ' '.join(map(str,cutout))
						curr_max_line = line_count
				else:
					g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
					h.write(str(baseline) + ' ' + curr_max_cutout + '\n')
					print_count += 1
					curr_max_ch = ch
					curr_max_frame = frame
					curr_max_amp = amp
					curr_max_cutout = ' '.join(map(str,cutout))
					curr_max_line = line_count

			line_count += 1

		# include append for last line
		g.write(str(curr_max_ch) + ' ' + str(curr_max_frame) + ' ' + str(curr_max_amp) + '\n')
		h.write(str(baseline) + ' ' + curr_max_cutout + '\n')
		print_count += 1

	print("Number of filtered spikes =", print_count)



def amps(cutout, baseline):
	peak = np.min(np.asarray(cutout[6:16]))
	bl = np.median(np.asarray(cutout[:5]))
	amp = baseline - peak
	if amp > 0:
		return amp
	elif bl - peak > 0:
		return bl - peak
	else:
		return 0


def localisation(file, chpos, clen=26, medians=False):
	chpos = np.load(chpos)

	# coord arrays
	X = []
	Y = []

	g = open('local_spikes', 'w')
	h = open('pca_cutouts', 'w')

	with open(file) as f:
		line_count = 0
		count2 = 0
		count_emptystring = 0
		for line in f:
			line_count += 1
			line = map(float, line.split())
			baseline = line[0]
			cutouts = line[1:]
			no_ch = int(len(cutouts))/clen
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
				cutout = np.array(cutouts[i*clen+1:(i+1)*clen])
				cutout_amp = amps(cutout, baseline)
				if chID == 384:
					print('384 line =', line_count)
				if cutout_amp > max_amp:
					max_amp = cutout_amp
					max_amp_cutout = ' '.join(map(str,cutout))
				amps_.append(cutout_amp)
				chpos_.append(chpos[chID])
			# Write max cutout to file
			if max_amp_cutout == '':
				# print "max_amp_cutout == '':" + str(amps_)
				count_emptystring += 1
				continue
			h.write(max_amp_cutout + '\n')
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
			else:
				print("Broken line =", line_count)
				print(cutouts)
				print(median, amps_)

		print("count2 =", count2)
		print("empty string =", count_emptystring)
