#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Thu 18. Dec 2025
#
# Purpose:      Controls the simulation by Euler-Maruyama temporal integration. 
# =============================================================================

import numpy as np
from models import cartpole_dynamics

def simulate(initial_state, controller, params,
             T=10.0, dt=0.01, noise_std=0.0):

    n_steps = int(T / dt)
    state = initial_state.copy()

    history = []

    for i in range(n_steps):
        t = i * dt

        theta = state[2]
        u = controller.compute(-theta, dt)

        dx = cartpole_dynamics(state, t, u, params, noise_std)
        state = state + dx * dt   # Euler

        history.append(state.copy())

    return np.array(history)

def main() -> None:
    pass

if __name__ == "__main__":
    main()

