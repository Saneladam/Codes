#!/usr/bin/env python3

# %% changes

# %% imports
import subprocess
import readline
import os

from ollama import chat

from rich.console import console
from rich.live import live
from rich.text import text
from dataclasses import dataclass, field

from commands import commands
from ollama_chat import *
from session import session, ask

console = console()


# %% function definition
commands = {}


def command(name):

    def wrapper(func):
        commands[name] = func
        return func

    return wrapper


def print_help():
    console.print(f"\t[model  : '{session.model}']")
    console.print(f"\t[i-file : '{session.ifile}']")
    console.print(f"\t[o-file : '{session.ofile}']")
    console.print()
    console.print("\t/q            quit")
    console.print("\t/h            help")
    # console.print("\t/c            clear conversation")
    # console.print("\t/m            change model")
    console.print("\t/i <file>     set input file")
    console.print("\t/o <file>     set output file")
    # console.print("\t/l            load log")
    # console.print("\t/s            edit system prompt")
    console.print("\t/y            copy last response")
    console.print("\t/! <command>  execute shell command")
    console.print()


def copy_function(text):
    subprocess.run(
        ["wl-copy"],
        input=text,
        text=true,
    )


def handle_command(session, command: str) -> bool:
    """
    returns true if command was handled.
    """
    if command == "/q":
        raise systemexit
    if command == "/h":
        print_help()
        return true
    # if command == "/c":
    #     messages[:] = [messages[0]]
    #     console.print("conversation cleared.")
    #     return true
    # if command == "/l":
    #     console.print("load conversation (todo)")
    #     return true
    if command == "/y":
        try:
            copy_function(session.last_answer)
            console.print("[green]last message copied into clipboard.[/]")
        except nameerror as e:
            console.print(f"[red]upsi, no message to be copied. :/[/red]")
        except exception as e:
            console.print(f"[red]error at copy function: {e}[/red]")

        return true
    if command.startswith("/i "):
        session.ifile = command[3:].strip()
        console.print("input file:", session.ifile)
        return true
    if command.startswith("/o "):
        session.ofile = command[3:].strip()
        console.print("output file:", session.ofile)
        return true
    if command.startswith("/!"):
        subprocess.run(
            command[2:],
            shell=true,
            executable=os.environ["shell"],
        )
        return true

    return false


def ask(session, prompt: str):
    session.messages.append(
        {
            "role": "user",
            "content": prompt,
        }
    )
    answer = ""
    try:
        stream = chat(model=session.model, messages=session.messages, stream=true)
        text = text()
        with live(text, refresh_per_second=20, console=console) as live:
            for chunk in stream:
                piece = chunk["message"]["content"]
                answer += piece
                text.append(piece)
                live.update(text)
    except exception as e:
        console.print(f"[red]{e}[/red]")
        return
    session.messages.append(
        {
            "role": "assistant",
            "content": answer,
        }
    )
    session.last_answer = answer


# %% main
def main():
    # console.rule("[bold red] L O C A L I A")
    print("-" * 33)
    print_help()
    print("-" * 33)
    while True:
        try:
            prompt = console.input("[bold blue]>>> [/]")
        except (EOFError, KeyboardInterrupt):
            console.print("\n[red]Interrupted.[/]")
            break
        if not prompt:
            continue
        if prompt.startswith("/"):
            handle_command(session, prompt)
            continue
        ask(session, prompt)


if __name__ == "__main__":
    main()
