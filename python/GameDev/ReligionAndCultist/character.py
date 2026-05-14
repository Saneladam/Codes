#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 01. Feb 2026
#
# Purpose:      Manage the behaviour of the character 
# =============================================================================

import pygame

class Character:
    """A class to manage the character."""
    def __init__(self, ai_game):
        self.screen = ai_game.screen
        self.screen_rect = ai_game.screen.get_rect()

        # Load the character image and get its rect.
        self.image = pygame.image.load('images/priest.png')
        self.rect = self.image.get_rect()

        # Start each new character at the bottom center of the screen.
        self.rect.midbottom  = self.screen_rect.midbottom

    def blitme(self):
        """Draw the character at tits current location.""" 
        self.screen.blit(self.image, self.rect) 
