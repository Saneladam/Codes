#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 08. Mar 2026
#
# Purpose:      This agent will tailor higly specific prompts for Cabo and save
#               its output.
# =============================================================================

from ollama import generate
import subprocess
from pathlib import Path
import sys

OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

ALFEREZ_PROMPT = """
You are ALFEREZ.

Mission:
Create a precise execution plan prompt for CABO, an agentic code and text wirter that will execute the orders prompted to generate a file.

CABO can only receive one instruction at a time and produces one output file.

Rules:

- Break the objective into independent tasks
- Each task must produce ONE output file
- Specify TEXT or CODE based on the requirements for the plan.
- Instructions must be precise and work as independent block codes
- No explanations No markdown No tables No headers
- Output ONLY task lines
- Generate a modest number of tasks
- Check that the plan follows the format strictly. If not, regenerate the plan.
- The tasks must tackle the plan in a methodic and good-practices and higly specific set of instructions that will be transmitted to CABO
- Your role is to act as a good Task Manager and provide tailored instructions to CABO

STRICT FORMAT (one task per line):

TASK01|TEXT|filename.txt|instruction
TASK02|CODE|script.py|instruction
TASK03|TEXT|notes.txt|instruction
...

Requirements:

- No spaces around the |
- No additional text before or after
- Tasks must be sequential

When generating CODE tasks:
- The instruction must explicitly say:
  "Write a complete runnable Python script"
  or
  "Write a Python function"
- The instruction must specify:
  - function names
  - inputs
  - outputs
  - libraries allowed
- The instruction must end with:
  "Output only the code."
"""

MODEL = "llama3.2:latest"


def create_plan(objective):
    prompt = f"""
{ALFEREZ_PROMPT}
Objective:
{objective}
"""
    response = generate(model=MODEL, prompt=prompt, options={"temperature": 0.2})
    return response["response"]


def parse_tasks(plan):
    tasks = []
    for line in plan.splitlines():
        if not line.startswith("TASK"):
            continue
        parts = line.split("|")
        task_id = parts[0].strip()
        mode = parts[1].strip()
        output_name = parts[2].strip()
        instruction = parts[3].strip()
        tasks.append(
            {
                "task": task_id,
                "mode": mode,
                "output_name": output_name,
                "instruction": instruction,
            }
        )
    return tasks


def call_cabo(instruction, mode):
    if mode == "CODE":
        flag = "-C"
    else:
        flag = "-T"
    process = subprocess.run(
        ["python", "Cabo.py", flag], input=instruction, text=True, capture_output=True
    )
    return process.stdout


def save_output(name, content):
    path = OUTPUT_DIR / name
    with open(path, "w") as f:
        f.write(content)


def execute_mission(objective):
    plan = create_plan(objective)
    print("####################################################################")
    print("ALFEREZ PLAN:")
    print(plan)
    print("####################################################################")
    save_output("ALFEREZ_PLAN.txt", plan)
    tasks = parse_tasks(plan)
    for task in tasks:
        print(f"Executing {task['task']} {task['output_name']}")
        instruction = f"""
MISSION CONTEXT:
{objective}

TASK:
{task['instruction']}
"""
        result = call_cabo(instruction, task["mode"])
        # Guardar output de Cabo con task_id
        output_filename = f"{task['output_name']}"
        save_output(output_filename, result)

        # Guardar la instrucción original
        instruction_filename = f"{task['task']}_INSTRUCTION.txt"
        save_output(instruction_filename, task["instruction"])


def main():
    if len(sys.argv) < 2:
        print("Usage: python Alferez.py 'mission'")
        sys.exit()
    objective = sys.argv[1]
    execute_mission(objective)

    mission_file = OUTPUT_DIR / "MISSION.txt"
    mission_file.write_text(objective)
    print("####################################################################")
    print("Launching CAPORAL review")
    subprocess.run(["python", "Caporal.py"])
    print("####################################################################")


if __name__ == "__main__":
    main()
