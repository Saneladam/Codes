#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Mon 29. Jun 2026
#
# Purpose:      Chatting with the local AI
# =============================================================================

from ollama import chat
from ollama import ChatResponse

response: ChatResponse = chat(
    model="xentriom/gemma-4-12B-coder-fable5-composer2.5-v1",
    messages=[
        {
            "role": "user",
            "content": "You are a proficient coder in python."
            + "I want you to write a ollama test generator with model:'xentriom/gemma-4-12B-coder-fable5-composer2.5-v1' that is fluent in python and bash so that I can input a script and it will enhance the code with comments, more pythonic lenguage and more robust implementation and comments on the right spots, or alternatively I would enter an input that depicts what i want the AI to build, lastlty it has to be able to support both modes in combination, so that i can give the code and also an instruction."
            + "It is important to only respond with the curated code, ready to be executed in python.",
        },
    ],
)
print(response["message"]["content"])
# or access fields directly from the response object
# print(response.message.content)
