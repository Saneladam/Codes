#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Use removesuffix() to clean and get the name of the .txt.
#               In this case .md
# =============================================================================

import os 

def main() -> None:
    print("On my Notes directory I have:")
    print("-----------------------------")
    for note in os.listdir("/home/akash/Documents/Notes"):
        print(f"\t{note.removesuffix('.md')}")
    print("-----------------------------")

if __name__ == "__main__":
    main()
