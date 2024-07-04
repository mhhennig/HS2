from pathlib import Path
from typing import TypedDict, Union


class Params(TypedDict):
    bandpass: bool
    freq_min: float
    freq_max: float
    chunk_size: int
    rescale: bool
    rescale_value: float
    common_reference: str
    spike_duration: float
    amp_avg_duration: float
    threshold: float
    min_avg_amp: float
    AHP_thr: float
    neighbor_radius: float
    inner_radius: float
    peak_jitter: float
    rise_duration: float
    decay_filtering: bool
    decay_ratio: float
    localize: bool
    save_shape: bool
    out_file: Union[str, Path]
    left_cutout_time: float
    right_cutout_time: float
    verbose: bool


DEFAULT_PARAMS: Params = {
    'bandpass': True,
    'freq_min': 300.0,
    'freq_max': 6000.0,

    'chunk_size': 100000,

    'rescale': True,
    'rescale_value': -1280.0,

    'common_reference': 'average',

    'spike_duration': 1.0,
    'amp_avg_duration': 0.4,
    'threshold': 10.0,
    'min_avg_amp': 5.0,
    'AHP_thr': 0.0,

    'neighbor_radius': 90.0,
    'inner_radius': 70.0,

    'peak_jitter': 0.2,
    'rise_duration': 0.26,

    'decay_filtering': False,
    'decay_ratio': 1.0,

    'localize': True,

    'save_shape': True,
    'out_file': 'HS2_detected',
    'left_cutout_time': 0.3,
    'right_cutout_time': 1.8,

    'verbose': True
}
