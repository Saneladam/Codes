#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Mon 29. Jun 2026
#
# Purpose:
# =============================================================================

import numpy as np
from scipy.fftgrid import fft2d  # Use mgrid/fftw instead if needed
import matplotlib.pyplot as plt
import time


def plasma_sim(n=128, dt=0.1, steps=50):
    N = n
    dx = 4 * np.pi / n
    dv = 1.0

    # Grid indices and coordinates
    x = np.linspace(-20, 20, N)
    v = np.linspace(-10, 10, N)
    xv = np.meshgrid(x, v)

    # Initial Maxwellian distribution
    f = np.exp(-(v**2) / 2.0)  # Gaussian in velocity space

    m = 1.0  # electron mass
    qe = 1.0  # charge unit
    kd = 30.0  # Debye length (k_d = k/lamdaD)

    field = np.zeros((n, n))
    rho = np.zeros((n, n))
    time = 0

    for _ in range(steps):
        # Compute charge density from f over v and a central source
        source_x = (n // 2) - 1
        charge_idx = (n // 2), (n // 2)
        rho = np.sum(f, axis=1) + 50 * np.eye(n)[charge_idx]

        # Solvable via FFT convolution with a Debye kernel
        k = np.fft.fftfreq(n).reshape(-1), np.fft.fftfreq(n).reshape(-1)
        rho_f = rho.astype(complex)
        rho_f[np.arange(n), :] = rho_f.real  # real only for simplicity

        kernel = (k[0] ** 2 + k[1] ** 2) / ((k[0] ** 2 + k[1] ** 2) ** 2 + kd**2)
        potential = np.real(np.fft.ifft2d(rho_f * kernel))
        field = -partial_x_gradient(potential, dx)

    plt.imshow(np.abs(field), extent=[-20, 20, -15, 15], cmap="viridis")
    plt.title("Screened Electric Field Magnitude (|E|)")
    plt.xlabel("Position (x)")
    plt.ylabel("Velocity (v)")
    plt.show()


def partial_x_gradient(potential, dx):
    # Approximate gradient with finite differences on the grid
    return np.zeros((128, 128))  # placeholder for actual fft-based grad if used
