"""Approximate solid-curve readings from Heggie--Aarseth Fig.4, not author data."""
import numpy as np

# Scan-reading uncertainty roughly 5--10%; these points are not exact targets.
REFERENCE = np.array([[1, 3.6, .6], [3, 5.43, 2.0], [10, 10.8, 9.5],
                      [30, 3.13, 3.53], [100, .35, .39]])


def sample(data, times):
    if min(times)<data[0, 1] or max(times)>data[-1, 1]:
        raise ValueError("history does not cover requested times; no extrapolation")
    return np.column_stack([np.interp(times, data[:, 1], data[:, c]) for c in (2, 3)])
