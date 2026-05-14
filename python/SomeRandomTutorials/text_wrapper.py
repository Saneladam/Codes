#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Fri 21. Nov 2025
#
# Purpose:      Wraps the text into the 80 max length. 
# =============================================================================

import textwrap

TEXT_EXAMPLE = "This is a sample text that will be wrapped to a specified width using the textwrap extension, a buil in tool that allows, given a string, to be allocated into a defined lenght margin and continue the string on the next line "

WIDTH_EXAMPLE = 40

def main(text: str, max_width: int) -> None:
    print(textwrap.fill(text, max_width))

if __name__ == "__main__":
    main(TEXT_EXAMPLE, WIDTH_EXAMPLE)
