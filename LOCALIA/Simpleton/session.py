#!/usr/bin/env python3

# %% Imports
import subprocess
import readline
import os

from ollama import chat

from dataclasses import dataclass, field

MODEL = "xentriom/gemma-4-12B-coder-fable5-composer2.5-v1"

SYSTEM_PROMPT = (
    "You are a machine with precise experty in software engineer, Linux,"
    "scientific high performance programmer and technical writer."
    "Produce concise, correct and practical computer output written in bash or python3 "
    "with focus on giving an efficient and sensible response that can be executed "
    "by a python or bash interpreter. Always include the proper shebang when generating scripts."
)


@dataclass
class Session:

    model: str = MODEL

    ifile: str = ""
    ofile: str = ""

    last_answer: str = ""

    messages: list = field(
        default_factory=lambda: [{"role": "system", "content": SYSTEM_PROMPT}]
    )


session = Session()


def ask(session, prompt: str):
    session.messages.append(
        {
            "role": "user",
            "content": prompt,
        }
    )
    answer = ""
    try:
        stream = chat(model=session.model, messages=session.messages, stream=True)
        text = Text()
        with Live(text, refresh_per_second=20, console=console) as live:
            for chunk in stream:
                piece = chunk["message"]["content"]
                answer += piece
                text.append(piece)
                live.update(text)
    except Exception as e:
        console.print(f"[red]{e}[/red]")
        return
    session.messages.append(
        {
            "role": "assistant",
            "content": answer,
        }
    )
    session.last_answer = answer
