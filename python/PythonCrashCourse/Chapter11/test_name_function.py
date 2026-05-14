#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:      Test for name function
# =============================================================================

import pytest

from name_function import get_formatted_name


def test_first_last_name() -> None:
    """Do names like Janis Joplin work?"""
    formatted_name = get_formatted_name("janis", "joplin")
    assert formatted_name == "Janis Joplin"


def test_first_middle_last_name() -> None:
    """Do names like Wolfgang Amadeus Mozart work?"""
    formatted_name = get_formatted_name("wolfgang", "mozart", "amadeus")
    assert formatted_name == "Wolfgang Amadeus Mozart"
