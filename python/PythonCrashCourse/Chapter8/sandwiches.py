#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Arguments and functions
# =============================================================================

from profile import build_profile

def sandwich(*toppings):
    print("Your sandwich has the following ingredients:")
    for topping in toppings:
        print(f"\t- {topping}")
    print("That's all.")

def main() -> None:
    user = build_profile("Román", "García", "Guill", allergies="none", field="nuclear")
    print(f"Hello {user['first_name'].title()}")
    print(user)
    sandwich("mozzarella", "mustard", "ham", "albahaca")

    user = build_profile("Alberto", "Muñoz", allergies="cheese", field="paramedic")
    print(f"Hello {user['first_name'].title()}")
    print(user)
    sandwich("onion", "ham", "tomato")
    pass


if __name__ == "__main__":
    main()
