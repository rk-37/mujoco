"""Interceptor control strategies.

Each module in this package exposes a single function:

    control(frame: np.ndarray, sensors: dict) -> np.ndarray  # (4,) motor thrusts

Select one at run time:  run_intercept.py --controller <module_name>
(e.g. --controller hover). Drop a new file here to add a new strategy.
"""
