#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Try it for your self.
# =============================================================================


def RentalCar():
    kind_of_car = input("What kind of car do you want?\n>")
    print(f"\nLet me see if I can find something like a {kind_of_car}.")


def RestaurantSeating():
    num_people = int(input("How many are you?\n>"))
    if num_people >= 8:
        print("\nIt looks like you will need to wait a little bit ...")
    else:
        print(f"\nThis way please, table for {num_people}.")


def MultipleOfTen():
    divisor = int(input("Write a number to see if it is divisible by 10\n"))
    if divisor % 10 == 0:
        print("Yes it was")
    else:
        print("No, baka")
    name = input("What is your name")
    if divisor % 10 == 0:
        print("Yes it was")
    else:
        print("No, baka")


def Sandwich():
    orders = {
        "Natalie": "double-cheese",
        "Clara": "single-cheese",
        "Elena": "veggie",
        "Maria": "double-cheese",
        "Sandra": "jamón",
    }
    completed = []
    to_deliver = [value for value in orders.values()]
    while to_deliver:
        print(f"Still need to deliver: \n{to_deliver}\n")
        just_coocked = input("What is that chef?\n> ").strip()
        new_plate = True
        for customer in [name for name in orders.keys()]:
            if (orders[customer] == just_coocked) and (new_plate):
                to_deliver.remove(just_coocked)
                completed.append(customer)
                print(
                    f"\nDingDingDing! {customer}!\nYour {just_coocked} sandwich is ready"
                )
                new_plate = False
        if new_plate:
            print("\nSorry chef you mistook the order")
        elif not to_deliver:
            print("\nThat was all! Lets go for a smoke")
        else:
            print("\nGreat! Lets continue")


def main() -> None:
    # RentalCar()
    # RestaurantSeating()
    # MultipleOfTen()
    Sandwich()
    pass


if __name__ == "__main__":
    main()
