# Information on the clustering methods in HS2

The clustering interface in HS2 is quite generic, so different algorithms can through a common set of commands. The code will combine spike locations (x,y) with a feature vector (any dimension), and cluster this combined data set.

## General workflow

We assume we have a `herdingspikes` instance with detected spikes ready to go, and this object is called ``H``. To see how to get this, see the [Information om detection](detection.md).

First we create an instance of the ``HSClustering`` class, passing ``H``, and storing the object reference in ``C``:

```python
from hs2 import HSClustering
C = HSClustering(H)
```

Next, we can call a function to compute features. Currently available are PCA

```python
C.ShapePCA(pca_ncomponents=2, pca_whiten=True);
```

and ICA

```python
C.ShapeICA(ica_ncomponents=2, ica_whiten=True);
```

The dimensionality is set with the parameter ``p/ica_ncomponents``.

Now, we can simply call the clustering function with appropriate parameters:

```python
H.CombinedClustering(alpha=alpha, clustering_algorithm=clustering_algorithm, cluster_subset=cluster_subset)
```

This has three important parameters:
* ``alpha`` is the scaling of the features, relative to the locations. This has to be adjusted to the dimensions used for  locations, which may be channel numbers or metric distances.
* ``clustering_algorithm`` points to the sklearn.cluster class to be used, and defaults to sklearn.cluster.MeanShift.
* ``cluster_subset`` is the number of data points actually used to form clusters. In a large recording (say 10s of millions of spikes), it is very inefficient to cluster every single event. Instead, the clusters can be formed using only ``cluster_subset`` events, and then all events are assigned to their nearest cluster centre. A good choice is around ``cluster_subset=100000``.

Further parameters are passed to the clustering function.

Finally, the result can be saved into a hdf5 container:

```python
C.SaveHDF5(["sorted.hdf5"])
```

## Using Mean Shift

This is the default clustering method in HS2, and is based on the implementation available in [Scikit-learn](http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html) (parallel code contributed by Martino Sorbaro). The version included [here](../clustering) is modified to guarantee good performance for large data sets on multiple CPU cores.

The command to use this is:

```python
C.CombinedClustering(cluster_subset=cluster_subset, alpha=alpha, bandwidth=bandwidth, bin_seeding=False, min_bin_freq=1, n_jobs=-1)
```

Specific parameters:
*  `bandwidth` - Kernel size.
*  `bin_seeding` - If `True`, the data is binned before clustering - this gives substantial speed improvements, but can make the result worse. If `cluster_subset` is used to limit the number of data points used to form clusters, this can be set to `False`.
*  `min_bin_freq` - Minimum number of spikes per bin (only if `bin_seeding` is `True`).
*  `n_jobs` - Number of threads to run simultaneously (for all available cores, use `-1`).

## Using (H)DBSCAN

Two other density based methods are [DBSCAN](http://scikit-learn.org/stable/modules/generated/sklearn.cluster.dbscan.html) and [HDBSCAN](https://hdbscan.readthedocs.io/en/latest/comparing_clustering_algorithms.html). These give good results, and can unlike Mean Shift, detect and remove noise events. HDBSCAN appears particularly suited, and scales well. To use this, install the `hdbscan` package with `pip` or `conda`.

Example usage:

```python
from sklearn.cluster import DBSCAN
C.CombinedClustering(alpha=0.4, clustering_algorithm=DBSCAN,
                    min_samples=5, eps=0.2, n_jobs=-1)
```
