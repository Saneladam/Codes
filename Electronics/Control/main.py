#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Thu 18. Dec 2025
#
# Purpose:      Ensables the inverse pendulum. 
# =============================================================================

import matplotlib
matplotlib.use("Agg")

import numpy as np
from controllers import PIDController
from simulation import simulate
from visualization import animate_cartpole

params = {
    "m": 0.1,
    "M": 1.0,
    "l": 1.0,
    "g": 9.81
}

initial_state = np.array([0.0, 0.0, 0.2, 0.0])

controller = PIDController(kp=50, ki=0.0, kd=10)

history = simulate(
    initial_state,
    controller,
    params,
    T=10.0,
    dt=0.001,
    noise_std=0.5
)

animate_cartpole(
    history,
    l=params["l"],
    step=10
)

def main() -> None:
    pass

if __name__ == "__main__":
    main()

