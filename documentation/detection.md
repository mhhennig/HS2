# Information on the HS2 spike detection method

Spike detection is based on a very fast, online-capable method developed by Oliver Muthmann described
[here](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Following detection, spatial spike locations are estimated, and redundant events removed. The results are stored in a flat binary file for fast access.

---

Assuming we have a created a ``Probe`` object called ``P``, use:

```python
H = HSDetection(Probe, to_localize, cutout_start, cutout_end, threshold,
                maa=0, maxsl=12, minsl=3, ahpthr=0, out_file_name=out_file, save_all=False)
```

Parameters:
* to_localize - If set to `False`, spikes will only be detected, not localised. Note the data cannot be sorted then.
* cutout_start - Number of frames to save backwards from spike peak.
* cutout_end - Number of frames to save forward from spike peak.
*  threshold - Detection threshold, this is given in multiples of the estimated noise level. Note the variability measure is not a variance, but a percentile much better suited for the highly non-Gaussian noise in the data. As a rule of thumb, this can be set to around 10.
*  maa - Minimal average amplitude of the event peak.
*  maxsl - Dead time in frames after spike peak, used for further testing.
*  minsl - Number of frames used for determining average spike amplitude.
*  ahpthr - The signal should go below this value within maxsl-minsl frames
* out_file_name - Base file name (without extension) for the output file(s).
* save_all - If `True`, store additional debug output.

---

Now we can read the result like so:

```python
H.LoadDetected()
```
