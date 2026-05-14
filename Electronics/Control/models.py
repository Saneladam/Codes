#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Thu 18. Dec 2025
#
# Purpose:      Creates the dinamics fo the pendulum 
# =============================================================================

import numpy as np

def cartpole_dynamics(state, t, u, params, noise_std=0.0):
    """
    state = [x, x_dot, theta, theta_dot]
    u     = fuerza aplicada al carro
    """

    x, x_dot, theta, theta_dot = state

    m = params["m"]      # masa del péndulo
    M = params["M"]      # masa del carro
    l = params["l"]      # longitud del péndulo
    g = params["g"]

    # Ruido estocástico (fuerza externa tipo turbulencia)
    noise = noise_std * np.random.randn()

    total_force = u + noise

    sin_t = np.sin(theta)
    cos_t = np.cos(theta)

    denom = M + m * sin_t**2

    x_ddot = (total_force + m * sin_t * (l * theta_dot**2 + g * cos_t)) / denom
    theta_ddot = (-total_force * cos_t - m * l * theta_dot**2 * cos_t * sin_t
                   - (M + m) * g * sin_t) / (l * denom)

    return np.array([x_dot, x_ddot, theta_dot, theta_ddot])



def main() -> None:
    pass

if __name__ == "__main__":
    main()

