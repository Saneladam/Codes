#!/usr/bin/env python3

# %% CHANGES
"""
- [ ] Make the /! commands execute in a specific shell command
- [ ] Make the handeling of the Errno 111 Connection refused and an integration for the ollama run service
- [ ] Make a input and output command
- [ ] Have a list of possibilities for different Model, and agents.
"""

# %% Imports
import subprocess
import readline
import os

from ollama import chat

from rich.console import Console
from rich.live import Live
from rich.text import Text
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

    messages: list = field(default_factory=lambda: [
        {
            "role": "system",
            "content": SYSTEM_PROMPT
        }
    ])
session = Session()

# %% Function Definition
COMMANDS = {}

def command(name):

    def wrapper(func):
        COMMANDS[name] = func
        return func

    return wrapper

def print_help():
    console.print(f"\t[Model  : '{session.model}']")
    console.print(f"\t[I-File : '{session.ifile}']")
    console.print(f"\t[O-File : '{session.ofile}']")
    console.print()
    console.print("\t/q            Quit")
    console.print("\t/h            Help")
    # console.print("\t/c            Clear conversation")
    # console.print("\t/m            Change model")
    console.print("\t/i <file>     Set Input file")
    console.print("\t/o <file>     Set Output file")
    # console.print("\t/l            Load log")
    # console.print("\t/s            Edit system prompt")
    console.print("\t/y            Copy last response")
    console.print("\t/! <command>  Execute shell command")
    console.print()


def copy_function(text):
    subprocess.run(
        ["wl-copy"],
        input=text,
        text=True,
    )


def handle_command(session, command: str) -> bool:
    """
    Returns True if command was handled.
    """
    if command == "/q":
        raise SystemExit
    if command == "/h":
        print_help()
        return True
    # if command == "/c":
    #     messages[:] = [messages[0]]
    #     console.print("Conversation cleared.")
    #     return True
    # if command == "/l":
    #     console.print("Load conversation (TODO)")
    #     return True
    if command == "/y":
        try:
            copy_function(session.last_answer)
            console.print("[green]Last message copied into clipboard.[/]")
        except NameError as e:
            console.print(f"[red]Upsi, no message to be copied. :/[/red]")
        except Exception as e:
            console.print(f"[red]Error at copy function: {e}[/red]")

        return True
    if command.startswith("/i "):
        session.ifile = command[3:].strip()
        console.print("Input file:", session.ifile)
        return True
    if command.startswith("/o "):
        session.ofile = command[3:].strip()
        console.print("Output file:", session.ofile)
        return True
    if command.startswith("/!"):
        subprocess.run(
            command[2:],
            shell=True,
            executable=os.environ["SHELL"],
        )
        return True

    return False


def ask(session, prompt: str):
    session.messages.append(
        {
            "role": "user",
            "content": prompt,
        }
    )
    answer = ""
    try:
        stream = chat(
            model=session.model,
            messages=session.messages,
            stream=True
        )
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


# %% Main
def main():
    # console.rule("[bold red] L O C A L I A")
    print("-"*33)
    print_help()
    print("-"*33)
    while True:
        try:
            prompt = console.input("[bold blue]>>> [/]")
        except (EOFError, KeyboardInterrupt):
            console.print("\n[red]Interrupted.[/]")
            break
        if not prompt:
            continue
        if prompt.startswith("/"):
            handle_command(session,prompt)
            continue
        ask(session,prompt)


if __name__ == "__main__":
    main()
