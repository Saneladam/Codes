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
lA = 1  # m
lB = 1  # m
mA = 1  # kg located at the end point of the pendulum
mB = 1  # kg located at the end point of the pendulum

timeFINAL = 50  # s
thetaA0 = 0.90 * np.pi
thetaB0 = 0.80 * np.pi
FPS = 60
frame_step = int((1 / FPS) / dt)

# %% Initial Values
times = np.arange(0, timeFINAL, dt)  # time of the simulation
xAs = np.zeros_like(times)
yAs = np.zeros_like(times)
xBs = np.zeros_like(times)
yBs = np.zeros_like(times)
thetaAs = np.zeros_like(times)
thetaBs = np.zeros_like(times)
omegaAs = np.zeros_like(times)
omegaBs = np.zeros_like(times)

# %% Iterable
for i in range(len(times) - 1):
    """
    Cambiar ecuaciones de doble pendulo
    """
    omegaAs[i] = (
        omegaAs[i - 1] + (g / lA) * np.sin(thetaAs[i - 1]) * dt if i != 0 else 0
    )
    omegaBs[i] = (
        omegaBs[i - 1] + (g / lB) * np.sin(thetaBs[i - 1]) * dt if i != 0 else 0
    )
    thetaAs[i] = thetaAs[i - 1] + omegaAs[i] * dt if i != 0 else thetaA0
    thetaBs[i] = thetaBs[i - 1] + omegaBs[i] * dt if i != 0 else thetaB0

    xAs[i] = lA * np.sin(thetaAs[i])
    yAs[i] = -lA * np.cos(thetaAs[i])
    xBs[i] = lB * np.sin(thetaBs[i]) + xAs[i]
    yBs[i] = -lB * np.cos(thetaBs[i]) + yAs[i]


# %% Table registration
path = OUTPUT_DIR / "double_pendulum.dat"
with open(path, "w") as f:
    f.write("# times , xAs , yAs , xBs , yBs , omegaAs , omegaBs , thetaAs , thetaBs ")
    for i in range(len(times)):
        f.write(
            " "
            + str(times[i])
            + " "
            + str(xAs[i])
            + " "
            + str(yAs[i])
            + " "
            + str(xBs[i])
            + " "
            + str(yBs[i])
            + " "
            + str(omegaAs[i])
            + " "
            + str(omegaBs[i])
            + " "
            + str(thetaAs[i])
            + " "
            + str(thetaBs[i])
        )
        f.write("\n")

# %% Plotting
fig, ax = plt.subplots()
ax.set_xlim(-1.2 * lA - 1.2 * lB, 1.2 * lA + 1.2 * lB)
ax.set_ylim(-1.2 * lB - 1.2 * lB, 1.2 * lA + 1.2 * lB)
ax.set_aspect("equal")
ax.set_title("Simple pendulum (Euler integration)")
time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

(lineA,) = ax.plot([], [], "o-", color="blue", lw=2)
(lineB,) = ax.plot([], [], "o-", color="red", lw=2)


def update(i):
    xA = [0, xAs[i]]
    yA = [0, yAs[i]]
    xB = [xAs[i], xBs[i]]
    yB = [yAs[i], yBs[i]]
    t = times[i]
    time_text.set_text(f"t = {t:.2f} s")
    lineA.set_data(xA, yA)
    lineB.set_data(xB, yB)
    return (lineA, lineB, time_text)


frames = range(0, len(times), frame_step)

ani = FuncAnimation(fig, update, frames=frames, interval=1000 / FPS, blit=True)
plt.show()

# %% End message
print("All Done")
