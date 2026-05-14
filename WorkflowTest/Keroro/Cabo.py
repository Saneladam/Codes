#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sat 07. Mar 2026
#
# Purpose:      Follows instructions as precise as possible.
# =============================================================================

import sys
import argparse

from ollama import generate

MODEL_TEXT = "llama3.2:latest"
MODEL_CODE = "codellama:latest"

CABO_TEXT_PROMPT = """
You are CABO-T.

Mission:
Execute the user instruction with maximum precision.

Rules:

- Be concise
- Be exact
- No unnecessary explanation
- No preamble
- No commentary
- No formatting unless required
- Output only the result

Priority order:

1 Precision
2 Correctness
3 Clarity
4 Brevity
"""
CABO_CODE_PROMPT = """
You are CABO-C.

Mission:
Write correct, efficient, and elegant code.

Rules:
- Output ONLY code
- No explanation
- No markdown
- No comments unless necessary
- Prefer simple and efficient solutions
- Code must run

Priority order:

1 Correctness
2 Efficiency
3 Simplicity
4 Elegance
"""


def parse_args():
    parser = argparse.ArgumentParser(description="CABO AI executor")
    parser.add_argument(
        "-C",
        "--code",
        action="store_true",
        help="Code generation mode",
    )
    parser.add_argument(
        "-T",
        "--text",
        action="store_true",
        help="Text generation mode (default)",
    )
    return parser.parse_args()


def run_llm(prompt: str, model: str):
    response = generate(
        model=model,
        prompt=prompt,
        options={
            "temperature": 0.1,
            "top_p": 0.9,
        },
    )
    return response["response"]


def main() -> None:
    args = parse_args()
    task = sys.stdin.read().strip()
    if not task:
        sys.exit("No input provided.")
    if args.code:
        model = MODEL_CODE
        system_prompt = CABO_CODE_PROMPT
    else:
        model = MODEL_TEXT
        system_prompt = CABO_TEXT_PROMPT

    prompt = f"""
{system_prompt}

Task:
{task}

Output:
"""
    result = run_llm(prompt, model)
    print(result)


if __name__ == "__main__":
    main()
