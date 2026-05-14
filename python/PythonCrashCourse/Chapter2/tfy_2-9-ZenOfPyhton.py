#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Write addition, sutraction, multiplication and diviosn operation that result in the number 8.
# =============================================================================

import numpy as np

def main() -> None:
    EXTENT=20
    number_options = np.arange(0,EXTENT+1)
    print(number_options)
    FAVORITE_NUMBER = [8]
    for chosen_number in FAVORITE_NUMBER:
        print("=============================") 
        print(f"CHOSEN NUMBER: {chosen_number}")
        print("-----------------------------") 
        print("Checking sumation") 
        for a_sum in number_options:
            for b_sum in number_options:
                if a_sum + b_sum == chosen_number:
                    print(f"Can be wirtten as {a_sum} + {b_sum}.")
        print("-----------------------------") 
        print("Checking substraction")
        for a_sum in number_options:
            for b_sum in number_options:
                if a_sum - b_sum == chosen_number:
                    print(f"Can be wirtten as {a_sum} - {b_sum}.")
        print("-----------------------------") 
        print("Checking multiplication") 
        for a_sum in number_options:
            for b_sum in number_options:
                if a_sum * b_sum == chosen_number:
                    print(f"Can be wirtten as {a_sum} * {b_sum}.")
        print("-----------------------------")
        print("Checking division") 
        for a_sum in number_options:
            for b_sum in number_options:
                if b_sum == 0:
                    continue
                elif a_sum / b_sum == chosen_number:
                    print(f"Can be wirtten as {a_sum} / {b_sum}.")
        print("=============================")

if __name__ == "__main__":
    main()
