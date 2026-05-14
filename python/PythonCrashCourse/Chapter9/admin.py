#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 25. Jan 2026
#
# Purpose:      Exercises using classes
# =============================================================================

from user import User

class Priviliges:
    def __init__(self, priviliges):
        self.priviliges = priviliges

    def show_power(self):
        print("As an administrator you can:")
        for privilige in self.priviliges:
            print(f"\t-{privilige}")

class Administrator(User):
    def __init__(
        self,
        username,
        password,
        priviliges=Priviliges(["Be God", "Own everything"]),
    ) -> None:
        super().__init__(username, password, login_attempts=1)
        username = self.name.upper()
        self.priviliges = priviliges
