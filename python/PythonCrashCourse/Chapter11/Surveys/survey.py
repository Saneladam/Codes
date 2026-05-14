#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 27. Jan 2026
#
# Purpose:       
# =============================================================================

class AnonymousSurvey:

    def __init__(self,question) -> None:
        self.question = question
        self.responses=[]

    def show_question(self):
        print(self.question)

    def store_response(self,new_response):
        self.responses.append(new_response)
    
    def show_result(self):
        print("Survey results:")
        for response in self.responses:
            print(f"-\t{response}")
