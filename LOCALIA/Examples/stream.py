#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Mon 29. Jun 2026
#
# Purpose:      Chatting with the local AI
# =============================================================================

from ollama import chat

# message_content = "You are a proficient coder. I want you to write a makefile that runs a python script whenever there is a change inside the Data directory."
message_content = "Define in the context of plasma Physics what are the ways to create new REs, such as avalanche, Compton scattering, tritium decay, hot-tail, dreicer."

stream = chat(
    model="gemma3",
    messages=[{"role": "user", "content": message_content}],
    stream=True,
)

for chunk in stream:
    print(chunk["message"]["content"], end="", flush=True)
print("This is the response message content")
