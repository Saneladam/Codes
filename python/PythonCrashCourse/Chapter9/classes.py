#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Exercises using classes
# =============================================================================


def doggie():
    class Dog:
        """An attempt to model a dog."""

        def __init__(self, name, age):
            """Constructor."""
            self.name = name
            self.age = age

        def sit(self):
            """Simulate a dog sitting in response to a command."""
            print(f"{self.name} is now sitting.")

        def roll_over(self):
            """Simulate a dog rolling over in response to a command."""
            print(f"{self.name} rolled over!")

    # def do_tricks(tricks):
    #     while True:
    # while True:
    #     trick = input("\t> ").strip().lower()
    #     if trick in tricks.keys():

    my_dog = Dog("Yeti", 17)

    print(f"My dog's name is {my_dog.name}.")
    print("It can do some tricks like: sit")
    my_dog.sit()
    print("It can do some tricks like: roll over")
    my_dog.roll_over()


def icecreamstands():

    from restaurant import Restaurant, IceCreamStand

    restaurant1 = Restaurant("Hermanos Pollo", "Chicken rotiserie")
    restaurant1.describe_restaurant()
    restaurant2 = Restaurant("Capitán Almejas", "seafood")
    restaurant2.describe_restaurant()
    restaurant3 = Restaurant("Tres Copas", '"tapas"')
    restaurant3.describe_restaurant()
    icecreamstand = IceCreamStand(
        "Primos Frigo", ["choco-mint", "cacahuete", "dulce de leche"]
    )
    icecreamstand.describe_restaurant()
    icecreamstand.describe_flavours()
    icecreamstand = IceCreamStand("Manolo")
    icecreamstand.describe_flavours()


def admin():

    from user import User
    from admin import Administrator

    first_user = User("Roman", "abc123")
    first_user.greet_user()
    first_user.ask_password()
    first_user = Administrator(first_user.name, first_user.password)
    print(f"{first_user.name} is now an Administrator\n")
    first_user.priviliges.show_power()


def electric_car():

    from car import Car, ElectricCar

    my_coal = Car("Xara Picaso", 2004)
    my_coal.get_specs()
    print()
    my_leaf = ElectricCar("XiYinPing", 3024)
    my_leaf.get_specs()
    print("------------------")
    my_leaf.battery.upgrade_battery()
    my_leaf.get_specs()


def dice():

    import random

    class Dice:
        """Descripción de la clase"""

        def __init__(self, sides=6):
            """Constructor."""
            self.sides = sides

        def roll_die(self):
            roll = random.randint(1, self.sides)
            print(roll)

    all_sides = [6, 10, 20]
    for case_side in all_sides:
        case_dice = Dice(case_side)
        print(f"\nRolling dice with {case_side} sides")
        i = 0
        while i < 6:
            i += 1
            case_dice.roll_die()


def lottery():

    import random

    destiny = (
        "N",
        "S",
        "V",
        "X",
        "Y",
        "P",
        "R",
        "L",
        "D",
        "M",
        "6",
        "4",
        "9",
        "7",
        "2",
    )

    selection = []
    my_ticket = []
    number_of_games = 0
    while len(my_ticket) < 4:
        my_ticket.append(random.choice(destiny))
    print(f"My lucky ticket is {my_ticket}")
    while selection != my_ticket:
        number_of_games += 1
        selection = []
        while len(selection) < 4:
            selection.append(random.choice(destiny))
    print(f"I just had to play {number_of_games} times. :))")


def main() -> None:
    # icecreamstands()
    # admin()
    # electric_car()
    # dice()
    lottery()
    pass


if __name__ == "__main__":
    main()
