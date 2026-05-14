#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Thu 18. Dec 2025
#
# Purpose:      Creates the matplotlib visualization 
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os

def animate_cartpole(history, l=1.0, outdir="frames", step=10):
    os.makedirs(outdir, exist_ok=True)

    x = history[:, 0]
    theta = history[:, 2]

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.set_xlim(-2, 2)
    ax.set_ylim(-1.5, 1.5)

    for i in range(0, len(x), step):

        if i % (10 * step) == 0:
            print(f"Rendering frame {i}/{len(x)}")

        ax.cla()
        ax.set_xlim(-2, 2)
        ax.set_ylim(-1.5, 1.5)

        ax.plot([x[i]], [0], "ks", markersize=8)

        px = x[i] + l * np.sin(theta[i])
        py = l * np.cos(theta[i])
        ax.plot([x[i], px], [0, py], "r-", lw=2)

        ax.axis("off")

        plt.savefig(f"{outdir}/frame_{i:05d}.png", dpi=100)

    plt.close(fig)

def main() -> None:
    pass

if __name__ == "__main__":
    main()

