#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Wed 23. Sep 2026
#
# Purpose:      Builds a 3d shape builder.
# =============================================================================

import numpy as np

# =============================================================================
# Base object
# =============================================================================


def sphere(radius=1.0, n_phi=200, n_theta=400):
    phi = np.linspace(0, np.pi, n_phi)
    theta = np.linspace(0, 2 * np.pi, n_theta)

    phi, theta = np.meshgrid(phi, theta, indexing="ij")

    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)

    return x, y, z, theta, phi


# =============================================================================
# Radial transformations
# =============================================================================


def radial_transform(object_3d, delta_r):
    """
    Apply a radial deformation to the object.

    delta_r is an array containing the radial displacement for
    every point of the surface.
    """

    x, y, z, theta, phi = object_3d

    radius = np.sqrt(x**2 + y**2 + z**2)

    new_radius = radius + delta_r

    # Avoid division by zero.
    direction_x = x / radius
    direction_y = y / radius
    direction_z = z / radius

    return (
        direction_x * new_radius,
        direction_y * new_radius,
        direction_z * new_radius,
        theta,
        phi,
    )


# =============================================================================
# 1. Waves
# =============================================================================


def wave(object_3d, amplitude=0.2, frequency=8.0):
    """
    Create a sinusoidal radial wave around the object.
    """

    x, y, z, theta, phi = object_3d

    delta_r = amplitude * np.sin(frequency * theta) * np.sin(phi)

    return radial_transform(object_3d, delta_r)


# =============================================================================
# 2. Spikes
# =============================================================================


def spikes(object_3d, amplitude=0.5, frequency=12.0, power=4.0):
    """
    Create radial spikes.

    power controls how sharp the spikes are.
    """

    x, y, z, theta, phi = object_3d

    pattern = np.abs(np.sin(frequency * theta) * np.sin(frequency * phi))

    pattern = pattern**power

    delta_r = amplitude * pattern

    return radial_transform(object_3d, delta_r)


# =============================================================================
# 3. Dents
# =============================================================================


def dents(
    object_3d,
    amplitude=0.5,
    theta0=0.0,
    phi0=np.pi / 2,
    sigma=0.25,
):
    """
    Create a local depression on the surface.

    theta0, phi0 define the position of the dent.
    sigma controls its size.
    """

    x, y, z, theta, phi = object_3d

    # Angular distance from the centre of the dent.
    dtheta = np.angle(np.exp(1j * (theta - theta0)))

    dphi = phi - phi0

    distance_squared = dtheta**2 + dphi**2

    delta_r = -amplitude * np.exp(-distance_squared / (2 * sigma**2))

    return radial_transform(object_3d, delta_r)


# =============================================================================
# Object construction
# =============================================================================

INITIAL_OBJECT = sphere()


OBJECT = INITIAL_OBJECT

OBJECT = wave(
    OBJECT,
    amplitude=0.15,
    frequency=8,
)

OBJECT = spikes(
    OBJECT,
    amplitude=0.35,
    frequency=10,
    power=5,
)

OBJECT = dents(
    OBJECT,
    amplitude=0.45,
    theta0=1.0,
    phi0=1.2,
    sigma=0.20,
)


FINAL_OBJECT = OBJECT


# =============================================================================
# Output
# =============================================================================

if __name__ == "__main__":
    x, y, z, theta, phi = FINAL_OBJECT

    np.savez(
        "object.npz",
        x=x,
        y=y,
        z=z,
        theta=theta,
        phi=phi,
    )
