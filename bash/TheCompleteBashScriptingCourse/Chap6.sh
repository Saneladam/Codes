#!/bin/bash

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Fri 26. Jun 2026
#
# Purpose:      Runs the Chapter 6 
# =============================================================================

first_name='roman'
last_name='garcia'

full_name="$first_name $last_name"

if [[ $first_name == 'roman' ]]; then
    echo "oh hey it's roman"
fi

echo "Hello $full_name!"
if ! (return 2>/dev/null); then
   : 
fi
