"""full_power -- example alternative controller.

Commands maximum thrust on all motors. With the drone level this just rockets it
upward -- useful only to confirm that swapping controllers changes the behavior.
"""

import numpy as np

CTRL_MAX = 6.0


def control(frame: np.ndarray, sensors: dict) -> np.ndarray:
    return np.full(4, CTRL_MAX, dtype=float)
