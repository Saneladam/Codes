#!/usr/bin/env python3

# %% Imports
import subprocess
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
    console.print(f"Model : {MODEL}")
    console.print()
    console.print("/q            Quit")
    console.print("/h            Help")
    console.print("/c            Clear conversation")
    console.print("/l            Load log")
    console.print("/y            Copy last response")
    console.print("/i <file>     Input file")
    console.print("/o <file>     Output file")
    console.print("/r <file>     Input + Output")
    console.print("/s            Edit system prompt")
    console.print("/m            Change model")
    console.print("/! <command>  Execute shell command")
    console.print()


def handle_command(command: str) -> bool:
    """
    Returns True if command was handled.
    """
    if command == "/q":
        raise SystemExit
    if command == "/h":
        print_help()
        return True
    if command == "/c":
        messages[:] = [messages[0]]
        console.print("Conversation cleared.")
        return True
    if command == "/l":
        console.print("Load conversation (TODO)")
        return True
    if command == "/y":
        console.print("Copy last answer (TODO)")
        return True
    if command.startswith("/i "):
        console.print("Input file:", command[3:])
        return True
    if command.startswith("/o "):
        console.print("Output file:", command[3:])
        return True
    if command.startswith("/r "):
        console.print("Input/Output:", command[3:])
        return True
    if command.startswith("/m "):
        global model
        model = command[3:].strip()
        console.print(f"Model: {model}")
        return True
    if command.startswith("/s "):
        console.print("System prompt updated.")
        return True
    if command.startswith("/!"):
        subprocess.run(command[2:].strip(), shell=True)
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
    except KeyboardInterrupt:
        console.print("\nInterrupted.")
        return
    except Exception as e:
        console.print(f"[red]{e}[/red]")
        return
    try:
        text = Text()
        with Live(text, refresh_per_second=20, console=console) as live:
            for chunk in stream:
                piece = chunk["message"]["content"]
                text.append(piece)
                live.update(text)
    except Exception as e:
        console.print(f"[red]{e}[/red]")
        return
    console.print()
    messages.append(
        {
            "role": "assistant",
            "content": answer,
        }
    )
    last_answer = answer


# %% Main
def main():
    console.rule("[bold red]Wellcome to LOCALIA")
    print_help()
    while True:
        try:
            console.print()
            prompt = console.input("\n[bold blue]>>> [/]")
        except (EOFError, KeyboardInterrupt):
            console.print()
            break
        if not prompt:
            continue
        if prompt.startswith("/"):
            handle_command(prompt)
            continue
        ask(prompt)


if __name__ == "__main__":
    main()
