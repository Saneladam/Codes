#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Exercises using classes
# =============================================================================

class User:
    def __init__(self, username, password, login_attempts=1) -> None:
        self.name = username
        self.password = password
        self.login_attempts = login_attempts

    def greet_user(self):
        print(f"\nWelcome {self.name} to the session.")

    def ask_password(self):
        while True:
            login_pass = input("Enter password:")
            if login_pass == self.password:
                print("Password correct\n")
                break
            elif self.login_attempts >= 3:
                print("Too many intents\n")
                break
            else:
                print("Incorrect, try again\n")
                self.login_attempts += 1
