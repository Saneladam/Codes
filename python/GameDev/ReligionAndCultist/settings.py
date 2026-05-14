#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 01. Feb 2026
#
# Purpose:      To store and manage the settings
# =============================================================================

class Settings:
    """A class to store all settings for Religion And Cultist."""
    
    def __init__(self):
        """Initialize the game's settings."""
        # Screen settings
        self.screen_width = 1200
        self.screen_height = 800
        #                 R   G   B
        # self.bg_color = (250,250,250)
        self.bg_color = (0,0,250)
