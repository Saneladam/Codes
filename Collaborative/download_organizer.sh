#!/bin/bash

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Fri 24. Jul 2026
#
# Purpose:      Solve Egor probelms. 
# =============================================================================

encouraging_short_sound() {
    pass
}

def_file_name=$( date ) 
week_name=$( date -V ) 
sound=encouraging_short_sound

dir_name=$( $week_name )
input_dir_name

if [[ -z $input_dir_name ]]; then
    pass
fi
mkdir -p "$week_name"
mkdir -p "$week_name/TEXT"
mkdir -p "$week_name/GEO"
find -n -i *.mtl -e { rm -- }
