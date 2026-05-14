#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 01. Feb 2026
#
# Purpose:      Manage the behaviour of the ship
# =============================================================================

import math

import pygame

class Ship:
    """A class to manage the ship."""

    def __init__(self, ai_game):
        self.screen = ai_game.screen
        self.settings = ai_game.settings
        self.screen_rect = ai_game.screen.get_rect()

        # Load the ship image and get its rect.
        self.image = pygame.image.load("images/ship.bmp")
        self.rect = self.image.get_rect()

        # Start each new ship at the bottom center of the screen.
        self.rect.midbottom = self.screen_rect.midbottom

        # Store a float for the ship's exact horizontal position.
        self.x = float(self.rect.x)
        self.y = float(self.rect.y)

        # Movement flag; start with a ship that's not moving
        self.moving_left = False
        self.moving_right = False
        self.moving_up = False
        self.moving_down = False

    def update(self):
        """Update the ship's position based on the movement flag."""
        if self.moving_left and self.rect.left > 0:
            if (self.moving_up or self.moving_down):
                self.x -= (self.settings.ship_speed)/math.sqrt(2)
            else:
                self.x -= self.settings.ship_speed
        if self.moving_right and self.rect.right < self.screen_rect.right:
            if (self.moving_up or self.moving_down):
                self.x += (self.settings.ship_speed)/math.sqrt(2)
            else:
                self.x += self.settings.ship_speed
        if self.moving_up and self.rect.top > 0:
            if (self.moving_right or self.moving_left):
                self.y -= (self.settings.ship_speed)/math.sqrt(2)
            else:
                self.y -= self.settings.ship_speed
        if self.moving_down and self.rect.bottom < self.screen_rect.bottom:
            if (self.moving_right or self.moving_left):
                self.y += (self.settings.ship_speed)/math.sqrt(2)
            else:
                self.y += self.settings.ship_speed

        self.rect.x = self.x
        self.rect.y = self.y

    def blitme(self):
        """Draw the ship at tits current location."""
        self.screen.blit(self.image, self.rect)
