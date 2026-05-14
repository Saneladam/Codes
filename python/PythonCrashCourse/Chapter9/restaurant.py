#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Exercises using classes
# =============================================================================


class Restaurant:
    def __init__(self, name, cuisine):
        self.name = name.title()
        self.type = cuisine.lower()

    def describe_restaurant(self):
        print(
            f"The renowed {self.name} restaurant is an exponenet of {self.type} cuisine."
        )


class IceCreamStand(Restaurant):
    def __init__(self, name, flavours=["vanilla", "strawberry", "chocolate"]):
        super().__init__(name, cuisine="gelateria")
        self.flavours = flavours

    def describe_flavours(self):
        print(f"{self.name} Ice Cream Stand has the following flavours:")
        for flavour in self.flavours:
            print(f"\t- {flavour}")
