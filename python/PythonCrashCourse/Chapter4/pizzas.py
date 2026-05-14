#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Pizzas, man i love them.
# =============================================================================

pizzas = ["Peperoni", "Quatri-Statzioni", "Calzone", "Salami", "Diavola", "Hawaian"]
my_favourite = pizzas[:]


def main() -> None:
    print("My favourite pizzas are:")
    print(my_favourite)
    # for fav in my_favourite:
    #     print(f"\t{fav}")
    print("Some people also have likings as:")
    pizzas.append("Caca")
    print(pizzas)

    print("\nBut I have to be honest,\nMy favourite pizzas are:")
    print(my_favourite)


if __name__ == "__main__":
    main()
