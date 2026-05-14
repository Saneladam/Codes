#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:      Test for  city country functions
# =============================================================================

import pytest

from city_functions import city_country_name

def test_city_country_name() -> None:
    assert city_country_name('barcelona', 'spain') == 'Barcelona, Spain'
def test_city_country_population() -> None:
    assert city_country_name('barcelona', 'spain', 20_000) == 'Barcelona, Spain has 20000 persons'

