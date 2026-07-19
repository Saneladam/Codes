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

from rich import prompt
from rich.console import Console
from rich.prompt import Prompt
from rich.live import Live
from rich.text import Text

console = Console()

# %% Constants

MODEL = "xentriom/gemma-4-12B-coder-fable5-composer2.5-v1"
model = MODEL
IFILE = ""
OFILE = ""

SYSTEM_PROMPT = (
    "You are an expert software engineer, Linux power user,"
    "scientific programmer and technical writer."
    "Poduce concise, correct and practical answers, with focus o "
    "pysics rigurosity, python and bash good code writing and technica "
    "poficiency in Linux open source KISS principle and"
)

messages = [
    {
        "role": "system",
        "content": SYSTEM_PROMPT,
    }
]


# %% Function Definition


def print_help():
    console.print(f"\t[Model  : '{MODEL}']")
    console.print(f"\t[I-File : '{IFILE}']")
    console.print(f"\t[O-File : '{OFILE}']")
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


def handle_command(command: str) -> bool:
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
            copy_function(last_answer)
            console.print("[green]Last message copied into clipboard.[/]")
        except NameError as e:
            console.print(f"[red]Upsi, no message to be copied. :/[/red]")
        except Exception as e:
            console.print(f"[red]Error at copy function: {e}[/red]")

        return True
    if command.startswith("/i "):
        global IFILE
        IFILE = command[3:].strip()
        # if IFILE:
        #     with open(IFILE,'write') as f:
        #
        console.print("Input file:", IFILE)
        return True
    if command.startswith("/o "):
        global OFILE
        OFILE = command[3:].strip()
        console.print("Output file:", OFILE)
        return True
    # if command.startswith("/m "):
    #     global model
    #     model = command[3:].strip()
    #     console.print(f"Model: {model}")
    #     return True
    # if command.startswith("/s "):
    #     console.print("System prompt updated.")
    #     return True
    if command.startswith("/!"):
        subprocess.run(
            command[2:],
            shell=True,
            executable=os.environ["SHELL"],
        )
        return True

    return False


def ask(prompt: str):
    global last_answer
    messages.append(
        {
            "role": "user",
            "content": prompt,
        }
    )
    answer = ""
    try:
        stream = chat(
            model=model,
            messages=messages,
            stream=True,
        )
        text = Text()
        with Live(text, refresh_per_second=20, console=console) as live:
            for chunk in stream:
                piece = chunk["message"]["content"]
                text.append(piece)
                live.update(text)
    except Exception as e:
        console.print(f"[red]{e}[/red]")
        return
    messages.append(
        {
            "role": "assistant",
            "content": answer,
        }
    )
    last_answer = answer


# %% Main
def main():
    console.rule("[bold magenta] L O C A L I A")
    print_help()
    while True:
        try:
            prompt = console.input("[bold blue]>>> [/]")
        except (EOFError, KeyboardInterrupt):
            console.print("\n[red]Interrupted.[/]")
            break
        if not prompt:
            continue
        if prompt.startswith("/"):
            handle_command(prompt)
            continue
        ask(prompt)


if __name__ == "__main__":
    main()
