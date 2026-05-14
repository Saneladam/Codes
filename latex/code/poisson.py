#!/usr/bin/env python3

# code/poisson.py
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 1, 200)
u = np.sin(np.pi * x)

np.save("data/u.npy", u)

plt.plot(x, u)
plt.xlabel("x")
plt.ylabel("u(x)")
plt.tight_layout()
plt.savefig("tex/figures/poisson.pdf")
