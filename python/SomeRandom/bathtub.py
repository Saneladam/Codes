#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Wed 19. Aug 2026
#
# Purpose:      Compare different numerical methods for solving the bathtub
#               filling problem.
# =============================================================================

cold_rate = 1.0  # bathtub / hour
hot_rate = 0.5  # bathtub / hour

desired_bathtubs = 1.0


def filling_rate():
    return cold_rate + hot_rate


def bathtub_volume(time):
    return filling_rate() * time


def residual(time):
    return bathtub_volume(time) - desired_bathtubs


# Bisection


def bisection(a, b, tolerance=1e-10, max_iterations=100):
    if residual(a) * residual(b) > 0:
        raise ValueError("The interval does not contain a root.")
    for i in range(max_iterations):
        midpoint = (a + b) / 2
        error = residual(midpoint)
        print(
            f"Iteration {i + 1:3d}: "
            f"t = {midpoint:.12f} h, "
            f"residual = {error:.3e}"
        )
        if abs(error) < tolerance:
            return midpoint
        if residual(a) * error < 0:
            b = midpoint
        else:
            a = midpoint
    return (a + b) / 2


# Newton-Raphson


def derivative(time):
    # dB/dt = total filling rate
    return filling_rate()


def newton(initial_guess, tolerance=1e-10, max_iterations=100):
    time = initial_guess
    for i in range(max_iterations):
        error = residual(time)
        print(
            f"Iteration {i + 1:3d}: " f"t = {time:.12f} h, " f"residual = {error:.3e}"
        )
        if abs(error) < tolerance:
            return time
        time = time - error / derivative(time)
    return time


# Main


def main():

    print("\n=== BISECTION ===\n")
    solution_bisection = bisection(a=0.0, b=1.0)
    print(f"\nSolution: {solution_bisection:.12f} h")
    print(f"         {solution_bisection * 60:.6f} min")
    print("\n=== NEWTON-RAPHSON ===\n")
    solution_newton = newton(initial_guess=1.0)
    print(f"\nSolution: {solution_newton:.12f} h")
    print(f"         {solution_newton * 60:.6f} min")


if __name__ == "__main__":
    main()
