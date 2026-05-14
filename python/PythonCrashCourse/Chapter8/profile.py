#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Arguments and functions
# =============================================================================

def build_profile(first, last, second_last=None, **user_info):
    user_info["first_name"] = first
    user_info["last_name"] = last
    if second_last:
        user_info["secondlast_name"] = second_last
    return user_info
