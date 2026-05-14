#!/usr/bin/env python3


# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:      Run thru the chapter
# =============================================================================

def get_formatted_name(first, second, middle=''):
    if middle:
        full_name = f"{first} {middle} {second}"
    else:
        full_name = f"{first} {second}"
    return full_name.title()
