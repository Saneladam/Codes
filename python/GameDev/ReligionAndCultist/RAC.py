#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 07. Feb 2026
#
# Purpose:      This is the main executer for the Religion and Cultist game. 
# =============================================================================

import sys

import pygame

from settings import Settings
from character import Character

class ReligionAndCultist:
    """Overall class to manage game assets and behaviour."""

    def __init__(self):
        """Initialize the game, and create game resources."""
        pygame.init()
        self.clock = pygame.time.Clock()
        self.settings = Settings()
        self.screen = pygame.display.set_mode((self.settings.screen_width,self.settings.screen_height))
        pygame.display.set_caption("Religion And Cultist")
        
        self.character = Character(self)

    def run_game(self):
        """Start the main loop for the game"""
        while True:
            # Watch for keyboard and mouse events.
            self._chech_events()
            self._update_screen()
            self.clock.tick(60)

    def _chech_events(self):
        """Respond to keypresses and mouse events."""
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()

    def _update_screen(self):
        """Update imgaes on the screen, and flip to new screen."""
        self.screen.fill(self.settings.bg_color)
        self.character.blitme()

        pygame.display.flip()


def main() -> None:
    rac = ReligionAndCultist()
    rac.run_game()

if __name__ == "__main__":
    main()

