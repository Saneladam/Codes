#!/bin/bash

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Mon 29. Jun 2026
#
# Purpose:      Test the trap command. 
# =============================================================================

cleanup() {
    echo "End function"
}

trap cleanup INT TERM EXIT

while true; do sleep 1; done
