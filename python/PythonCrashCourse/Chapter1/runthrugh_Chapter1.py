#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 24. Jan 2026
#
# Purpose:      Follow of Exampel for Chapter 2 hello  world 
# =============================================================================

date = "Today, loser"
first_name = "          román    "
first_name = first_name.rstrip()
first_name = first_name.lstrip()
last_name = "    garcía guill   "
last_name = last_name.rstrip()
last_name = last_name.lstrip()

full_name = f"{first_name.title()} {last_name.upper()}"

header_message = f"\t\t{date}\n"
hello_message = f"Hello, {first_name.title()}!\n"
goodbye_message = f"\nBye sucker!\nSigned:\t{full_name}\n"


def ask_body():
    body_message = input("Write the message here:\n")
    return body_message

def main() -> None:
    web_url="htps://roman_garcia_guill.com"
    print("-----------------------------")
    print("Welcome to Creating the email")
    content = ask_body()
    print("-----------------------------")
    print(f"{header_message}{hello_message} \n{content}\n{goodbye_message}")
    print("-----------------------------")
    print(f"Visit:\t{web_url.removeprefix('https://').removeprefix('htps://')}")

    
if __name__ == "__main__":
    main()
