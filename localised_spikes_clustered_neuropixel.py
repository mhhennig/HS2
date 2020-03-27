# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from IPython import get_ipython

# %% [markdown]
# # Download a data file from the UCL server
# 
# We're using a small [example data set](http://data.cortexlab.net/singlePhase3/) recorded by [Nick Steinmetz](http://www.nicksteinmetz.com/) at UCL.

# %%
#import urllib.request
#file_url = 'http://data.cortexlab.net/singlePhase3/data/rawDataSample.bin'
#file_name = "D:/Hopkins_20160722_g0_t0.imec.ap_CAR.bin"

#urllib.request.urlretrieve(file_url, file_name)


# %%
# Use of the `HSDetection` class


# %%
#import pyximport
#pyximport.install()
from herdingspikes.hs2 import HSDetection
from herdingspikes.probe import NeuroPixel
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# %%
# raw data file (as downloaded above)
data_path = "notebooks/rawDataSample.bin"

# detection parameters
file_directory = 'results/'
file_name = 'ProcessedSpikesNeuropixel'


# %%
Probe = NeuroPixel(data_path, neighbor_radius=120,
                   masked_channels=[36, 75, 112, 151, 188, 227, 264, 303, 340, 379, 384])

H = HSDetection(Probe, num_com_centers=2, threshold=20,
                left_cutout_time=0.4, right_cutout_time=1.0,
                maa=0, amp_evaluation_time=0.1, spk_evaluation_time=0.4,
                ahpthr=0, out_file_name=file_name, 
                file_directory_name=file_directory, decay_filtering=True, save_all=True)


# %%
Probe.show()


# %%
H.DetectFromRaw()


# %%
H.LoadDetected()


# %%
plt.figure(figsize=(10, 10))
H.PlotTracesChannels(15, window_size=100)


# %%
plt.figure(figsize=(14, 3))
ax = plt.subplot(111)
H.PlotAll(invert=True, s=1, alpha=0.2, ax=ax)
ax.set_xlim((0,600))

# %% [markdown]
# # Use of the `HSClustering` class

# %%
from herdingspikes.hs2 import HSClustering

# Load from file
#C = Clustering(['results/ProcessedSpikes_mea1k.bin', 'results/ProcessedSpikes_mea1k.bin'], cutout_length=41)

# Or if the spikes are already in memory
# simply load from the Detection class
C = HSClustering(H)


# %%
get_ipython().run_cell_magic('time', '', '# Compute features\nC.ShapePCA(pca_ncomponents=2, pca_whiten=True);')


# %%
#from sklearn.cluster import AgglomerativeClustering

#clusterer = AgglomerativeClustering(n_clusters=None, compute_full_tree=True, distance_threshold=50, linkage='single')
#from sklearn.cluster import DBSCAN
#dbscan = DBSCAN(min_samples=10, eps=1, n_jobs=-1, algorithm="ball_tree", leaf_size=30)
get_ipython().run_cell_magic('time', '', '# Cluster all spikes\n\n# This shows how a (potentially relatively small) subset of spikes \n# can be used to form clusters. All spikes are then assigned\n# to these clusters in batches. This allows clustering of very large \n# recordings (also multiple files) without exhausting memory.\n\nC.CombinedClustering(alpha=7, bandwidth=12.6, bin_seeding=False, n_jobs=-1, cluster_subset=100000)')
#get_ipython().run_cell_magic('time', '', '# Cluster all spikes\n\n# This shows how a (potentially relatively small) subset of spikes \n# can be used to form clusters. All spikes are then assigned\n# to these clusters in batches. This allows clustering of very large \n# recordings (also multiple files) without exhausting memory.\n\nC.CombinedClustering(alpha=0.8, clustering_algorithm=dbscan)')
#get_ipython().run_cell_magic('time', '', '# Cluster all spikes\n\n# This shows how a (potentially relatively small) subset of spikes \n# can be used to form clusters. All spikes are then assigned\n# to these clusters in batches. This allows clustering of very large \n# recordings (also multiple files) without exhausting memory.\n\nC.CombinedClustering(alpha=7, bin_seeding=False, n_jobs=-1, cluster_subset=100000, clustering_algorithm=clusterer)')


# %%
plt.figure(figsize=(20, 2))
ax = plt.subplot(111)
C.PlotAll(invert=True, s=1, ax=ax, max_show=100000, show_labels=False)
plt.title("MeanShift, bandwidth=5, no bin seeding")
ax.set_xlim((0,600))


# %%
units = range(16)
C.PlotShapes(units)


# %%
C.PlotNeighbourhood(9, radius=6, alpha=0.8)

# %%



# %%
# Save the results

C.SaveHDF5(file_name+"_sorted.hdf5")


# %%


