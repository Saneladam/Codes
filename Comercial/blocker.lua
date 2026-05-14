
-- =============================================================================
-- Authors:      Román García Guill
-- Contact:      romangarciaguill@gmail.com
-- Created:      Tue 31. Mar 2026
--
-- Purpose:      Manges the introduction of blocks in the nvim file.
--               Returns a block numbered in progression, if none then 0
-- =============================================================================

block="^# %%.*"
file=$1

grep "$block" "$file"
