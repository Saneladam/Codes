#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:      Run thru the chapter 
# =============================================================================
def example():
    def format_name(first,second):
        full_name = f"{first} {second}"
        return full_name.title()
    def get_names():
        print("Enter 'q' at any point to quit.\n")
        while True:
            first = input("What is your first name?\n\t")
            if first == 'q':
                break
            last = input("What is your last name?\n\t")
            if last == 'q':
                break
            formatted_name = format_name(first, last)
            print(f"\nYour name is \t{formatted_name}\n")
    get_names()


def main() -> None:
    example()
    pass

if __name__ == "__main__":
    main()

