#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 16. Jun 2026
#
# Purpose:      Generates a double pendulum.
#               First is just a normal pendulum
# =============================================================================

# %% Imports
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

OUTPUT_DIR = Path("./")

# %% Constants
dt = 0.001  # time definition in s
g = -9.81  # m/s**2
l = 1  # m
m = 1  # kg located at the end point of the pendulum

timeFINAL = 50  # s
theta0 = 0.50 * np.pi

FPS = 60
frame_step = int((1 / FPS) / dt)

# %% Initial Values
times = np.arange(0, timeFINAL, dt)  # time of the simulation
xs = np.zeros_like(times)
ys = np.zeros_like(times)
thetas = np.zeros_like(times)
omegas = np.zeros_like(times)

# %% Iterable
for i in range(len(times) - 1):
    omegas[i] = omegas[i - 1] + (g / l) * np.sin(thetas[i - 1]) * dt if i != 0 else 0
    thetas[i] = thetas[i - 1] + omegas[i] * dt if i != 0 else theta0

    xs[i] = l * np.sin(thetas[i])
    ys[i] = -l * np.cos(thetas[i])


# %% Table registration
path = OUTPUT_DIR / "simple_pendulum.dat"
with open(path, "w") as f:
    f.write("# times , xs , ys , omegas , thetas ")
    for i in range(len(times)):
        f.write(
            " "
            + str(times[i])
            + " "
            + str(xs[i])
            + " "
            + str(ys[i])
            + " "
            + str(omegas[i])
            + " "
            + str(thetas[i])
        )
        f.write("\n")

# %% Plotting
fig, ax = plt.subplots()
ax.set_xlim(-1.2 * l, 1.2 * l)
ax.set_ylim(-1.2 * l, 1.2 * l)
ax.set_aspect("equal")
ax.set_title("Simple pendulum (Euler integration)")
time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

(line,) = ax.plot([], [], "o-", lw=2)


def update(i):
    x = [0, xs[i]]
    y = [0, ys[i]]
    t = times[i]
    time_text.set_text(f"t = {t:.2f} s")
    line.set_data(x, y)
    return (line, time_text)


frames = range(0, len(times), frame_step)

ani = FuncAnimation(fig, update, frames=frames, interval=1000 / FPS, blit=True)
plt.show()

# %% Final message
print("All Done")
