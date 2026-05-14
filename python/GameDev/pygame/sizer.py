#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 21. Mar 2026
#
# Purpose:       
# =============================================================================

import pygame

pygame.init()
font = pygame.font.SysFont("IBM Plex Mono", 11)
text_surface = font.render("M", True, (0, 0, 0))
width, height = text_surface.get_size()
print(f"Ancho: {width}px, Alto: {height}px, Ratio: {height/width:.2f}")
