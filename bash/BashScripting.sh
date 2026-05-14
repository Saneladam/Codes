#!/bin/bash

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Thu 04. Dec 2025
#
# Purpose:      Go through the video: https://youtu.be/p0KKBmfiVl0 
#               Never say "If" writing a Bash script! (Exit codes & logical
#               operators) 
# =============================================================================
: << 'LOGICAL OPERATORS TO USE ELEGANTLY'
; && || &
LOGICAL OPERATORS TO USE ELEGANTLY

# 1) ;
echo "1) ;"
echo ""
echo "Hello there." ; echo "General Kenobi\!"
echo ""
echo "--------------------------------------------------------"
# 2) &&
# && is picky, only will run if the first command worked
echo "2) &&"
echo ""
cat .something && echo "This is something"
echo ""
echo "--------------------------------------------------------"
# 3) ||
# || will execute when the comand fails
echo "3) ||"
echo ""
cat .notathing || echo "This is not a thing"
echo ""
echo "--------------------------------------------------------"
# 4) &
# & will send the action to bg and execute the second one inmediately..
echo "4) &"
echo ""
sleep 5 & echo "Started the 5 seconds timer"
echo ""
echo "--------------------------------------------------------"

# =============================================================================

# condition1 && condition2 && condition3 || exit
