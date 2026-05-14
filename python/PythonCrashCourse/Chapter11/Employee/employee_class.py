#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Wed 28. Jan 2026
#
# Purpose:      Define the class
# =============================================================================


class Employee:
    """Descripción de la clase"""

    def __init__(self, first, last, salary):
        """Constructor."""
        self.first = first
        self.last = last
        self.salary = salary

    def give_rise(self, ammount=5000):
        self.salary += ammount

    def summary(self):
        print(f"{self.last.upper()}, {self.first.title()} earns {self.salary}")
