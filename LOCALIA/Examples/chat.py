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

message_content = "You are a proficient coder. I want you to write a makefile that runs a python script whenever there is a change inside the Data directory."

response: ChatResponse = chat(
    model="xentriom/gemma-4-12B-coder-fable5-composer2.5-v1",
    messages=[
        {
            "role": "user",
            "content": message_content,
        },
    ],
)
print("This is the response message content")

print(response["message"]["content"])
# or access fields directly from the response object
# print(response.message.content)
