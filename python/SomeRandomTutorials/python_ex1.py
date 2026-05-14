#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Fri 21. Nov 2025
#
# Purpose:      Follow along a python exercice focued on deepening python
#               understanding.
# =============================================================================

# Version 01: 
#       Dictionary that returns the letters of the words with 4 or more letters
def main() -> None:
    words = ["code", "python", "ai", "refactor", "bug"]
    length_map: dict[str, int] = {}
    for word in words:
        if len(word) > 4:
            length_map[word] = len(word)
    print(length_map)

# Version 02: 
#       Same but using comprehension
def main() -> None:
    words = ["code", "python", "ai", "refactor", "bug"]
    length_map = {word: len(word) for word in words if len(word) > 4}
    print(length_map)

# -----------------------------------------------------------------------------

# Version 03: 
#       Now also print the word with the index
def main() -> None:
    words = ["code", "python", "ai", "refactor", "bug"]
    for i in range(len(words)):
        print(i, words[i])

# Version 04: 
#       Same but in a more pythonic way (using python build in functions)
def main() -> None:
    words = ["code", "python", "ai", "refactor", "bug"]
    for i, word in enumerate(words):
        print(i, word)

# Version 05: 
#       Same but in a more pythonic way (using python build in functions)
def list_words(words: list[str]) -> None:
    for i, word in enumerate(words):
        print(i, word)
def main() -> None:
    words = ["code", "python", "ai", "refactor", "bug"]
    my_function = list_words
    my_function(words)

# Version 06: 
#       Same but in a more pythonic way (using python build in functions)
class ListWords:
    def __call__(self, words: list[str]) ->any:
        for i, word in enumerate(words):
            print(i, word)
def main() -> None:
    words = ["code", "python", "ai", "refactor", "bug"]
    list_words_fn = ListWords()
    list_words_fn(words)

# Version 07: 
#       Same but in a more pythonic way (using python build in functions)
class Transformer:
    def transform(self, my_list: list[str]) ->list[str]:
        return [word[::-1] for word in my_list] 
def do_something(transformer: Transformer, my_list: list[str]) -> None: 
    new_words = transformer.transform(my_list)
    print(new_words)
def main() -> None:
    words = ["code", "python", "ai", "refactor", "bug"]
    transformer = Transformer()
    do_something(transformer, words)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == "__main__":
    main()

