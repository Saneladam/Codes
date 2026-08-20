#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Tue 18. Aug 2026
#
# Purpose:
# =============================================================================

from ollama import chat

messages = [
    {
        "role": "user",
        "content": "Why is the sky blue?",
    },
]

for part in chat("gemma3", messages=messages, stream=True):
    print(part["message"]["content"], end="", flush=True)
