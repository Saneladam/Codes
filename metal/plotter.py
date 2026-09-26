#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Wed 23. Sep 2026
#
# Purpose:      Takes an image of the picture with a specific view point angle.
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

data = np.load("object.npz")

x = data["x"]
y = data["y"]
z = data["z"]


fig = plt.figure()

ax = fig.add_subplot(
    111,
    projection="3d",
)

ax.plot_surface(
    x,
    y,
    z,
    linewidth=0,
    antialiased=True,
)

ax.set_box_aspect((1, 1, 1))

plt.show()
