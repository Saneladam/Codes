#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Wed 28. Jan 2026
#
# Purpose:      Test for employee class
# =============================================================================

import pytest

from employee_class import Employee

@pytest.fixture
def jhon_doe():
    jhon = Employee('Jhon', 'Doe', 30000)
    return jhon

# @pytest.fixture
# def JuanGarcia():
#     juangarcia = Employee('Juan', 'García', '30000')
#     return juangarcia

def test_give_default_rise(jhon_doe):
    jhon_doe.give_rise()
    assert jhon_doe.salary == 35000

def test_give_custom_rise(jhon_doe):
    jhon_doe.give_rise(60)
    assert jhon_doe.salary == 30060
