"""hover -- placeholder controller.

Holds a constant hover thrust on all four motors. With the drone launched level,
this cancels gravity, so it flies a straight line in the throw direction and
misses. This is the baseline; real homing strategies go in sibling modules.
"""

import numpy as np

MASS = 0.272                 # kg (interceptor)
G = 9.81                     # m/s^2
HOVER_TOTAL = MASS * G       # ~2.67 N to hover
CTRL_MAX = 6.0               # per-motor thrust limit (model ctrlrange)


def control(frame: np.ndarray, sensors: dict) -> np.ndarray:
    """frame: (720,1280,3) uint8 RGB up-camera. sensors: {accel,gyro,alt}.
    Returns (4,) motor thrusts [px,nx,py,ny] in Newtons.

    TODO (homing): detect the target in `frame`, convert its bearing + IMU
    attitude into a desired lead direction (proportional navigation), and map
    that to differential per-motor thrust to re-aim and close. For now: hover.
    """
    cmd = np.full(4, HOVER_TOTAL / 4.0, dtype=float)
    return np.clip(cmd, 0.0, CTRL_MAX)
