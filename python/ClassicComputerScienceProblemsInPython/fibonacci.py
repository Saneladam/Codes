#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Wed 28. Jan 2026
#
# Purpose:      Fibonacci calculations
# =============================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# def fib1(n: int) -> int:
#     return fib1(n - 1) + fib1(n - 2)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# def fib2(n: int) -> int:
#     if n < 2:  # base case
#         return n
#     return fib2(n - 1) + fib2(n - 2)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# from typing import Dict

# memo: Dict[int, int] = {0: 0, 1: 1}


# def fib3(n: int) -> int:
#     if n not in memo:
#         memo[n] = fib3(n - 1) + fib3(n - 2)  # memoization
#         print(f"{n}\t->\t{memo[n]}")
#     return memo[n]


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# from functools import lru_cache


# @lru_cache(maxsize=None)
# def fib4(n: int) -> int:
#     if n < 2:  # base case
#         return n
#     return fib4(n - 1) + fib4(n - 2)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# def fib5(n: int) -> int:
#     if n == 0:
#         return n  # special case
#     last: int = 0  # initially set to fib(0)
#     next: int = 1  # initially set to fib(1)
#     for _ in range(1, n):
#         last, next = next, last + next
#     return next


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from typing import Generator


def fib6(n: int) -> Generator[int, None, None]:
    yield 0 # special case
    if n > 0: yield 1  # special case
    last: int = 0  # initially set to fib(0)
    next: int = 1  # initially set to fib(1)
    for _ in range(1, n):
        last, next = next, last + next
        yield next


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def main() -> None:
    for i in fib6(50):
        print(i)
    pass


if __name__ == "__main__":
    main()
