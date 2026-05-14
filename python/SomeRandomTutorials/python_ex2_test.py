#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Fri 21. Nov 2025
#
# Purpose:      To try to and add a simple test. 
# =============================================================================

import pytest

# Function we want to test
def is_prime(n: int) -> bool:
    if n <=1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def test_zero_is_prime() -> None:
    assert is_prime(0) == False
@pytest.mark.parametrize(
        "n, expected", [
            (2, True), 
            (3, True), 
            (4, False), 
            (5, True), 
            (6, False), 
            (7, True), 
            (8, False), 
            (17, True), 
            (0, False), 
            (-1, False)
        ])
def test_is_prime(n: int, expected: bool) -> None:
    assert is_prime(n) == expected
