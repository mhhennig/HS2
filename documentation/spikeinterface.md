# Using Herdingspikes with SpikeInterface

Herdingspikes is very easy to use with [SpikeInterface](https://spikeinterface.readthedocs.io/en/latest/). SpikeInterface provides a unified interface to spike sorting algorithms, and can be used to run Herdingspikes on data stored in a variety of formats.

First install the latest version of SpikeInterface from github:

```bash
git clone https://github.com/SpikeInterface/spikeinterface.git
cd spikeinterface
pip install -e .
```

Next, import SpikeInterface and create a `RecordingExtractor` object from the data. Here we use a `BioCamRecordingExtractor` that reads data from a BioCam recording. A list of supported formats is [here](https://spikeinterface.readthedocs.io/en/latest/modules/extractors.html#raw-data-formats), and documentation on the `RecordingExtractor` class is [here](https://spikeinterface.readthedocs.io/en/latest/modules/core.html#recordingextractor).

```python
import spikeinterface[full] as si

rx = si.BioCamRecordingExtractor('path/to/nwb/file.brw')
```

Pre-processing should not be required in many cases as Herdingspikes implements its own filters.
Now we can run the sorter:

```python
sx = si.run_sorter(
    "herdingspikes",
    rx,
    folder="hs_output",
    verbose=True,
    remove_existing_folder=True,
)
se.NpzSortingExtractor.write_sorting(sx, 'sorting.npz')
```

`sx` is a `SortingExtractor` object that contains the spike sorting results. The last line writes the sorting out to a file. Documentation on the `SortingExtractor` class is [here](https://spikeinterface.readthedocs.io/en/latest/modules/core.html#sorting).

If `remove_existing_folder=False`, the sorter will not overwrite the output folder if it already exists. If `verbose=True`, the sorter will print out more information about its progress.

## Setting sorter parameters

Parameters for the sorter can be set using a dictionary. A list of all parameters is [here](parameters.md#parameters].  Here we set the `clustering_min_bin_freq` parameter to 16:

```python
p = si.get_default_sorter_params("herdingspikes").copy()
p.update({"clustering_min_bin_freq":16})
sx = si.run_sorter(
    "herdingspikes",
    rx,
    folder="hs_output",
    verbose=True,
    remove_existing_folder=True,
    **p
)
```

## Drift correction with SpikeInterface

SpikeInterface provides sever drift correction algorithms that should be used for data from implanted probes. For a detailed documentation on drift correction, see [here](https://spikeinterface.readthedocs.io/en/latest/how_to/handle_drift.html).

To use drift correction with Herdingspikes, we simply apply a pre-processing step to the `RecordingExtractor` object before running the sorter:

```python

rx_drift_corrected, motion, motion_info = si.correct_motion(
    rx, preset="dredge", output_motion=True, output_motion_info=True
)
sx_drift_corrected = si.run_sorter(
    "herdingspikes",
    rx_drift_corrected,
    folder="hs_output_drift_corrected",
    verbose=True,
    remove_existing_folder=True,
)
```