#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Go through the execrices on Chapter 4 for loops and lists
# =============================================================================


def count_to_20(n=20) -> list:
    count_list = [number**3 for number in range(1, n + 1, 1)]
    return count_list


def main() -> None:
    count_to_onemillion = count_to_20(n=10)
    initial, final, sumation = (
        min(count_to_onemillion),
        max(count_to_onemillion),
        sum(count_to_onemillion),
    )
    print(f"\
The minimum is \t{initial}.\n\
The maximum is \t{final}.\n\
The sum is     \t{sumation}.")
    print(count_to_onemillion)

    print()
    print(f"The first 3 items in the list are: {count_to_onemillion[:3]}")
    print(f"The len of the list is {len(count_to_onemillion)}")
    print(
        f"The middle 3 items in the list are: {count_to_onemillion[int(len(count_to_onemillion)/2)-2:int(len(count_to_onemillion)/2)+1]}"
    )
    print(f"The last 3 items in the list are: {count_to_onemillion[-3:]}")


if __name__ == "__main__":
    main()
