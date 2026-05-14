#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Exercises using classes
# =============================================================================


class Car:
    def __init__(self, name, year) -> None:
        self.name = name
        self.year = year

    def get_specs(self):
        print(f"This car is a {self.name} from {self.year}.")


class Battery:
    def __init__(self, size) -> None:
        self.battery_size = size

    def upgrade_battery(self):
        if self.battery_size < 65:
            self.battery_size = 65
            print("Battery has been upgraded")

    def show_capacity(self):
        return self.battery_size

    def show_range(self):
        if self.battery_size > 60:
            self.range = self.battery_size * 28
        elif self.battery_size > 20:
            self.range = self.battery_size * 25
        else:
            self.range = self.battery_size * 22
        return self.range


class ElectricCar(Car):
    def __init__(self, name, year, battery_size=40) -> None:
        super().__init__(name, year)
        self.battery = Battery(battery_size)

    def get_specs(self):
        print(f"This electric car is a {self.name} from {self.year}.")
        print(
            f"The battery of {self.battery.show_capacity()} Ah can run up to {self.battery.show_range()} km."
        )
