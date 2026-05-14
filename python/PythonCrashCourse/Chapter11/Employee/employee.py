#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Wed 28. Jan 2026
#
# Purpose:      Playground for Employee
# =============================================================================

from employee_class import Employee

informatico = Employee("Juan", "Lorente", 26000)
informatico.summary()
print("Giving default rise")
informatico.give_rise()
informatico.summary()

informatica = Employee("Irene", "Lorente", 26000)
informatica.summary()
print("Giving default rise plus somthing else")
informatica.give_rise(7000)
informatica.summary()


def main() -> None:
    pass


if __name__ == "__main__":
    main()
