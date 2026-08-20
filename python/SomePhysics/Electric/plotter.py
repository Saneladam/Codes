#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
#
# Purpose:      Plot the electrostatic potential and electric field generated
#               by the potential_field.dat data file.
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

import scienceplots

# =============================================================================
# Configuration
# =============================================================================

DATA_FILE = "potential_field.dat"

NX = 101
NY = 101

# =============================================================================
# Load data
# =============================================================================

rows = []

with open(DATA_FILE, "r") as file:

    for line_number, line in enumerate(file, start=1):

        line = line.strip()

        # Ignore empty lines
        if not line:
            continue

        values = line.split()

        if len(values) == 3:

            x, y, potential_value = map(float, values)

            rows.append(
                [
                    x,
                    y,
                    potential_value,
                    np.nan,
                    np.nan,
                ]
            )

        elif len(values) == 5:

            x, y, potential_value, ex, ey = map(float, values)

            rows.append(
                [
                    x,
                    y,
                    potential_value,
                    ex,
                    ey,
                ]
            )

        else:

            raise ValueError(
                f"Unexpected number of columns at line " f"{line_number}: {len(values)}"
            )


data = np.array(rows)

x = data[:, 0].astype(int)
y = data[:, 1].astype(int)

V = data[:, 2]
Ex = data[:, 3]
Ey = data[:, 4]

# =============================================================================
# Convert flattened data into 2D arrays
# =============================================================================

potential = np.full((NX, NY), np.nan)

electric_x = np.full((NX, NY), np.nan)
electric_y = np.full((NX, NY), np.nan)


for i in range(len(data)):

    potential[x[i], y[i]] = V[i]

    if data.shape[1] == 5:
        electric_x[x[i], y[i]] = Ex[i]
        electric_y[x[i], y[i]] = Ey[i]


# Create mesh

X, Y = np.meshgrid(np.arange(NX), np.arange(NY), indexing="ij")


# =============================================================================
# Derived quantities
# =============================================================================

electric_field = np.sqrt(electric_x**2 + electric_y**2)


# =============================================================================
# Scientific plot style
# =============================================================================
plt.style.use(["science", "no-latex", "grid"])


# =============================================================================
# Figure 1: Electric potential
# =============================================================================
fig, ax = plt.subplots(figsize=(7, 5.5))
levels = np.linspace(np.nanmin(potential), np.nanmax(potential), 40)
contour = ax.contourf(X, Y, potential, levels=levels, cmap="coolwarm")

contour_lines = ax.contour(
    X,
    Y,
    potential,
    levels=levels[::4],
    linewidths=0.5,
)

ax.clabel(contour_lines, inline=True, fontsize=7, fmt="%.0f")

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_title(r"Electrostatic potential $V(x,y)$")

ax.set_aspect("equal")

cbar = fig.colorbar(contour, ax=ax)

cbar.set_label(r"$V$")

fig.tight_layout()

fig.savefig("potential.png", dpi=300, bbox_inches="tight")


# =============================================================================
# Figure 2: Electric field magnitude
# =============================================================================

fig, ax = plt.subplots(figsize=(7, 5.5))

levels = np.linspace(np.nanmin(electric_field), np.nanmax(electric_field), 40)

field_plot = ax.contourf(X, Y, electric_field, levels=levels, cmap="coolwarm")

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_title(r"Electric field magnitude $|\mathbf{E}|$")

ax.set_aspect("equal")

cbar = fig.colorbar(field_plot, ax=ax)

cbar.set_label(r"$|\mathbf{E}|$")

fig.tight_layout()

fig.savefig("electric_field_magnitude.png", dpi=300, bbox_inches="tight")


# =============================================================================
# Figure 3: Electric field lines
# =============================================================================

fig, ax = plt.subplots(figsize=(7, 5.5))

potential_plot = ax.contourf(X, Y, potential, levels=40, cmap="coolwarm")

ax.streamplot(
    X.T,
    Y.T,
    electric_x.T,
    electric_y.T,
    density=1.5,
    linewidth=0.8,
)

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_title(r"Electric potential and field lines")

ax.set_aspect("equal")

cbar = fig.colorbar(potential_plot, ax=ax)

cbar.set_label(r"$V$")

fig.tight_layout()

fig.savefig("potential_and_field.png", dpi=300, bbox_inches="tight")

plt.show()
