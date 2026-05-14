#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Mon 26. Jan 2026
#
# Purpose:      You already know what is this about
# =============================================================================

from pathlib import Path
import json


def pi_birthday():

    # path = Path('chapter_10/reading_from_a_file/pi_digits.txt')
    in_path = Path("chapter_10/reading_from_a_file/pi_million_digits.txt")
    out_path = Path("output.txt")
    try:
        content = in_path.read_text().rstrip()
    except FileNotFoundError:
        print("ERROR: missing input file")
        exit()
    pi_string = ""
    lines = content.splitlines()
    for line in lines:
        pi_string += line.strip()
    print(pi_string[:50])
    print(pi_string[-50:])
    print(len(pi_string))

    searcher = input(
        "Write the desired search number to look up at first million pi digits:\n>\t"
    )
    if searcher in pi_string:
        print("It is")
    else:
        print("It is not")

    out_path.write_text(searcher)


def Guest():
    guest_file = Path("guest.txt")
    assistants = ""
    while True:
        name = input("What is your name?\n(exit to quit)").strip().title()
        if name.lower() == "exit":
            break
        else:
            assistants += name + "\n"
    guest_file.write_text(assistants)


def division_calculator():
    print("Given two numbers a division of the first by the second will occur")
    print("Type 'q' to quit")
    while True:
        first_number = input("\nFirst number:\t")
        if first_number == "q":
            break
        second_number = input("Second number:\t")
        if second_number == "q":
            break
        try:
            answer = float(first_number) / float(second_number)
        except ZeroDivisionError:
            print("No division can be made with 0 as a denominator")
        except ValueError:
            print("You must enter numbers to do a division")
        else:
            print(answer)


def Cat_and_Dogs():
    cat_file = Path("cats.txt")
    dog_file = Path("dogs.txt")
    # cat_file = Path("cat.txt")
    # dog_file = Path("dog.txt")
    try:
        cat_content = cat_file.read_text()
    except FileNotFoundError:
        # print("Missing cat names :)")
        pass
    else:
        print("Some name of cats can be:")
        for name in cat_content.splitlines():
            print(f"- {name}")
    try:
        dog_content = dog_file.read_text()
    except FileNotFoundError:
        # print("Missing dog names :)")
        pass
    else:
        print("Some name of dogs can be:")
        for name in dog_content.splitlines():
            print(f"- {name}")


def CommonWords():
    book = Path("/home/akash/Lecture/Culture/BrothersKaramazov_FyodorDostoyevsky.txt")
    try:
        book_content = book.read_text().lower()
    except FileNotFoundError:
        print("Error on locating the book")
    else:
        unique_words = set(book_content.split())
        maximum_count = 0
        maximum_word = ""
        for word in unique_words:
            number_of_times = book_content.count(word + " ")
            print(f"{word}\t\t{number_of_times}")
            if number_of_times > maximum_count:
                maximum_count = number_of_times
                maximum_word = word
        print(
            f"\nTHE MOST COMMON WORD IS : {maximum_word.upper()} appearing {maximum_count} times."
        )


def number_writer():

    numbers = [x + 1 for x in range(10)]

    path = Path("numbers.json")
    content = json.dumps(numbers)
    path.write_text(content)


def number_reader():

    path = Path("numbers.json")
    content = path.read_text()
    numbers = json.loads(content)
    for number in numbers:
        print(f" -{number}")


def remember_me():
    path = Path("username.json")
    if path.exists():
        content = path.read_text()
        username = json.loads(content)
        print(f"Welcome back, {username}")
    else:
        username = input("What is your name?\n\t")
        content = json.dumps(username)
        path.write_text(content)
        print("\nWe will remember that for next time")


def FavouriteNumber():

    path = Path("fav_num_combined.json")
    if path.exists():
        content = path.read_text()
        number = json.loads(content)
        print(f"Your favourite number is {number}")
    else:
        number = input("what is your favourite number?\n")
        content = json.dumps(number)
        path.write_text(content)


def UserDictionary():

    def register_new_user_info(fields):
        user_dict={}
        for field in fields:
            user_dict[field]=input(f"What is your {field}?\n\t")
        content = json.dumps(user_dict)
        path.write_text(content)
        print("\nWe will remember that for next time")

    def get_users_info(user_dict):
        print("Welcome back")
        for field, respond in user_dict.items():
            print(f"I remember that your {field} was {respond}")

    def get_data(path):
        content = path.read_text()
        user_dict = json.loads(content)
        return user_dict

    def check_data(path):
        if path.exists():
            return True
        else:
            return False

    def check_for_user(name):
        confirmation = input(f"Are you {name}?\n")
        if confirmation == 'yes' or confirmation == 'y':
            return True
        else:
            return False


    path = Path("username_dictionary.json")
    fields = ["name", "age", "job"]

    if check_data(path):
        user_dict = get_data(path)
        name = user_dict['name']
        if check_for_user(name):
            get_users_info(user_dict)
        else:
            register_new_user_info(fields)

def main() -> None:
    # pi_birthday()
    # Guest()
    # division_calculator()
    # Cat_and_Dogs()
    # CommonWords()
    # number_writer()
    # number_reader()
    # remember_me()
    # FavouriteNumber()
    UserDictionary()
    pass


if __name__ == "__main__":
    main()
