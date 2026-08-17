#!/bin/bash

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Fri 24. Jul 2026
#
# Purpose:      Solve Egor probelms. 
# =============================================================================

dir_name=$( date +%V ) 

encouraging_short_sound() {
    pass
}

def_file_name=$( date +%d%m%y_%H%M ) 
input_dir_name() {
    read -r "input_dir"
    if [[ -z $input_dir ]]; then
        input_dir="$def_file_name"
    fi
}

input_dir_nameencouraging_short_sound
mkdir -p "$dir_name"
mkdir -p "$dir_name/TEXT"
mkdir -p "$dir_name/GEO"
find -n -i *.mtl -e { rm -- }
