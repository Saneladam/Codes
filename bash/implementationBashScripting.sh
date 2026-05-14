#!/bin/bash

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Thu 04. Dec 2025
#
# Purpose:      implement something on the video. 
# =============================================================================

# somehting done inelegantly
echo $EDITOR
if [ "$EDITOR" = "" ]; then
        EDITOR=vim
fi
echo $EDITOR

# now elegantly
echo $EDITOR
[ -z "$EDITOR" ] && EDITOR=vim
echo $EDITOR

