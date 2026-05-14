#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:      Run thru the chapter
# =============================================================================

from name_function import get_formatted_name

print("Enter 'q' at any point to quit.\n")
while True:
    first = input("What is your first name?\n\t")
    if first == "q":
        break
    middle = input("What is your middle name if any?\n\t")
    if middle == "q":
        break
    last = input("What is your last name?\n\t")
    if last == "q":
        break

    if not middle:
        formatted_name = get_formatted_name(first, last)
    else:
        formatted_name = get_formatted_name(first, last, middle)

    print(f"\nYour name is \t{formatted_name}\n")
