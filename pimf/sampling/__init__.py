"""
sampling doctring
"""
from .sampling import (
    IMFSample,
    IMFSampleList,
    draw_samples,
    smith2021_sampling,
    optimal_sampling
)

__all__ = ["IMFSample", "IMFSampleList", "draw_samples", "optimal_sampling", "smith2021_sampling"]