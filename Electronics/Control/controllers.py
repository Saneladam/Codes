#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Thu 18. Dec 2025
#
# Purpose:      Defines the classical PIDControllers 
# =============================================================================

class PIDController:
    def __init__(self, kp, ki, kd):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.integral = 0.0
        self.prev_error = 0.0

    def compute(self, error, dt):
        self.integral += error * dt
        derivative = (error - self.prev_error) / dt
        self.prev_error = error

        return (
            self.kp * error
            + self.ki * self.integral
            + self.kd * derivative
        )

def main() -> None:
    pass

if __name__ == "__main__":
    main()

