#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Use a variable to represent a person-s name, and then print that
#               person's name in lowercase, uppercase and title case.
# =============================================================================
name = "rodrigo diaz de vivar, EL CAMPEADOR"

def main() -> None:
    print(f"Input name:    \t{name}")
    print(f"Official name: \t{name.title()}")
    print(f"Uppercase name:\t{name.upper()}")
    print(f"Lowercase name:\t{name.lower()}")

if __name__ == "__main__":
    main()

