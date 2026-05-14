#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:       
# =============================================================================

from survey import AnonymousSurvey

question = "\nWhat language did you first learn to speak?"
language_survey = AnonymousSurvey(question)

language_survey.show_question()
print("Enter 'q' at any point to quit")
while True:
    response = input("Lenguage: ")
    if response == 'q':
        break
    language_survey.store_response(response)

print("\nThank you to everyone who participated!")
language_survey.show_result()
