#!/usr/bin/env python3

# =============================================================================
# Authors:      Roman Garcia Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 22. Nov 2025
#
# Purpose:      Explore the modules to use. 
# =============================================================================

"""
1)  Dataclasses
2)  Pathlib
3)  Functools
4)  Tomllib
5)  Graphlib raphlib
6)  Heapq
7)  Secrets
8)  Shutil
9)  Textwrap
10) Itertools
"""
# 1) dataclasses
from dataclasses import dataclass
@dataclass(frozen=True)
class Product:
    name: str
    price: float
    in_stock: bool = True
def main() -> None:
    product = Product(name="Widget", price=19.99)
    print(product)
# 2) pathlib
from pathlib import Path
def main() -> None:
    base = Path("my_project")
    config = base / "config" / "settings.toml"
    print("Config path:", config)
    if config.exists():
        print("File size:", config.stat().st_size)
    else:
        config.parent.mkdir(parents=True, exist_ok=True)
        config.write_text("[settings]\nname = 'Example'")
# 3) functools
from functools import cache, partial
@cache
def  power(base: int, exponent: int) -> int:
    print(f"Computing {base}^{exponent}")
    return base**exponent
def main() -> None:
    print("Calculating powers with caching:")
    print("Power of 2^10:", power(2,10))
    print("Power of 3^7:", power(3,7))
    print("Power of 4^5:", power(4,5))
    print("Power of 2^10 again (cached):", power(2,10))

    square = partial(power,exponent=2)
    cube = partial(power,exponent=3)

    print("Square of 5:", square(5))
    print("Square of 2:", cube(2))
    print("Square of 5 again (cached):", square(5))

# 4) tomllib
import tomllib
def main() -> None:
    with open("pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    print(data["project"]["description"])

# 5) graphlib
# 6) heapq
# 7) secrets
# 8) shutil
# 9) textwrap
#10) itertools


# ////////////////////////////////////////////////////////////////////////////
# def main() -> None:

if __name__ == "__main__":
    main()
