#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Find a quote from a famous person you admire,
#               print the quote and the name.
# =============================================================================

author = "Jiddu Krishnamurti"

quote = "Freedom from the desire for an answer is essential for the understanding of a problem"

def main() -> None:

    print(f'\n{author.title()} once said:\n\t"{quote}".')

if __name__ == "__main__":
    main()

