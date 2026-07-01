#!/bin/bash

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Fri 26. Jun 2026
#
# Purpose:       
# =============================================================================

my-func() {
    echo 'Hi'
    return 69
}

var=$(my-func)
code=$?

echo "output=$var, code=$code"
