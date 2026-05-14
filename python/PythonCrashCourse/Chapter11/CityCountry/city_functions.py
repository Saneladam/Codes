#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:      Some functions for the City Country code 
# =============================================================================

def city_country_name(city, country, population=None):
    if population:
        string =f"{city.title()}, {country.title()} has {population} persons"
    else: 
        string =f"{city.title()}, {country.title()}"

    return string
