#!/usr/bin/env python3

# =============================================================================
# Authors:      Román García Guill
# Contact:      romangarciaguill@gmail.com
# Created:      Sun 08. Mar 2026
#
# Purpose:      This will be Caporal, the validator and critic of the operation.
# =============================================================================

import subprocess
from pathlib import Path

from ollama import generate
import ast

MAX_ITERATIONS = 3

DANGEROUS_IMPORTS = {
    "os",
    "subprocess",
    "socket",
    "shutil",
}
DANGEROUS_FUNCTIONS = {
    "eval",
    "exec",
    "compile",
    "__import__",
}

OUTPUT_DIR = Path("outputs")

MODEL_TEXT = "llama3.2:latest"
MODEL_CODE = "codellama:latest"


CAPORAL_PROMPT_CODE = """
You are CAPORAL.

Mission:
Audit, validate and debug code produced by CABO.

Focus on:
- runtime errors
- security risks
- correctness
- good practices
- clarity
- optimization
- scientific rigor
- usefulness

Be strict but pragmatic.

Output format:

VERDICT: APPROVED or ITERATE

ANALYSIS:
technical explanation

IMPROVEMENTS:
specific fixes or NONE
"""


CAPORAL_PROMPT_TEXT = """
You are CAPORAL.

Mission:
Audit, validate and correct text produced by CABO.

Focus on:
- correctness
- clarity
- rigor
- usefulness
- completeness
- structure

Be strict but pragmatic.

Output format:

VERDICT: APPROVED or ITERATE

ANALYSIS:
technical explanation

IMPROVEMENTS:
specific fixes or NONE
"""


def security_scan(code):
    try:
        tree = ast.parse(code)
    except Exception as e:
        return False, f"Syntax error: {e}"
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for name in node.names:
                if name.name.split(".")[0] in DANGEROUS_IMPORTS:
                    return False, f"Dangerous import detected: {name.name}"
        if isinstance(node, ast.ImportFrom):
            if node.module and node.module.split(".")[0] in DANGEROUS_IMPORTS:
                return False, f"Dangerous import detected: {node.module}"
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                if node.func.id in DANGEROUS_FUNCTIONS:
                    return False, f"Dangerous function detected: {node.func.id}"
    return True, "OK"


def read_outputs():
    files = {}
    for path in OUTPUT_DIR.iterdir():
        if path.is_file():
            files[path.name] = path.read_text()
    return files


def read_mission():
    mission_file = OUTPUT_DIR / "ALFEREZ_PLAN.txt"
    if mission_file.exists():
        return mission_file.read_text()
    return ""

def read_task_instruction(task_id):
    task_file = OUTPUT_DIR / f"{task_id}_INSTRUCTION.txt"
    if task_file.exists():
        return task_file.read_text()
    return ""


def run_python(file_path):
    result = subprocess.run(
        ["python", "-I", "-B", str(file_path)],  # isolated mode  # no .pyc
        capture_output=True,
        text=True,
        timeout=5,
        cwd=file_path.parent,
        env={},
    )
    return {
        "stdout": result.stdout,
        "stderr": result.stderr,
        "returncode": result.returncode,
    }


def run_bash(file_path):
    result = subprocess.run(
        ["bash", str(file_path)], capture_output=True, text=True, timeout=10
    )
    return {
        "stdout": result.stdout,
        "stderr": result.stderr,
        "returncode": result.returncode,
    }


def evaluate_code(code, runtime, project_context):
    prompt = f"""
{CAPORAL_PROMPT_CODE}

PROJECT CONTEXT:
{project_context}

CODE UNDER REVIEW:
{code}

STDOUT:
{runtime["stdout"]}

STDERR:
{runtime["stderr"]}

Evaluate the code considering the whole project.
Detect duplicated logic or structural issues.
"""
    response = generate(model=MODEL_CODE, prompt=prompt, options={"temperature": 0.1})
    return response["response"]


def evaluate_text(text, project_context):

    prompt = f"""
{CAPORAL_PROMPT_TEXT}

PROJECT CONTEXT:
{project_context}

TEXT UNDER REVIEW:
{text}

Evaluate the text considering the whole project.
"""
    response = generate(model=MODEL_TEXT, prompt=prompt, options={"temperature": 0.1})
    return response["response"]


def process():
    mission = read_mission()
    files = read_outputs()
    for name, content in files.items():
        if name.endswith("_review.txt") or name == "ALFEREZ_PLAN.txt":
            continue
        path = OUTPUT_DIR / name

        # Extraer task_id del nombre del archivo
        # Espera formato: TASK_<task_id>_<original_filename>
        parts = name.split("_", 1)
        if len(parts) < 2:
            task_id = "UNKNOWN"
        else:
            task_id = parts[1]  # TASK01, TASK02, etc.
        task_instruction = read_task_instruction(task_id)

        iterations = 0
        while iterations < MAX_ITERATIONS:
            print(f"\nCAPORAL reviewing {name} (iteration {iterations+1})")
            project_context = build_project_context(read_outputs())
            # -------- CODE --------
            if name.endswith(".py"):
                safe, reason = security_scan(content)
                if not safe:
                    report = f"""
VERDICT: ITERATE
ANALYSIS:
Security violation detected.

DETAILS:
{reason}

IMPROVEMENTS:
Remove dangerous operations.
"""
                else:
                    runtime = run_python(path)
                    runtime["stdout"] = runtime["stdout"][:2000]
                    runtime["stderr"] = runtime["stderr"][:2000]
                    report = evaluate_code(content, runtime, project_context)
            # -------- BASH --------
            elif name.endswith(".sh"):
                report = evaluate_code(
                    content,
                    {"stdout": "", "stderr": "Execution skipped", "returncode": -1},
                    project_context,
                )
            # -------- TEXT --------
            else:
                report = evaluate_text(content, project_context)

            verdict = verdict_from_report(report)
            review_path = OUTPUT_DIR / f"{name}_review.txt"
            review_path.write_text(report)
            print(verdict)
            if verdict == "APPROVED":
                break
            print("Requesting fix from CABO")
            if name.endswith(".py") or name.endswith(".sh"):
                fixed = request_fix_code(content, report, mission, task_instruction)
            else:
                fixed = request_fix_text(content, report, mission, task_instruction)
            if not fixed.strip():
                print("CABO returned empty output")
                break
            path.write_text(fixed)
            content = fixed
            iterations += 1
        if iterations >= MAX_ITERATIONS:
            print(f"Max iterations reached for {name}")


def global_project_review(project_context):
    prompt = f"""
You are CAPORAL.

Evaluate the entire project.

PROJECT STRUCTURE:
{project_context}

Check for:

- duplicated logic
- missing files
- inconsistent design
- poor architecture

Output format:

PROJECT_VERDICT: APPROVED or ITERATE

ANALYSIS:
...
"""
    response = generate(model=MODEL_TEXT, prompt=prompt, options={"temperature": 0.1})
    return response["response"]


def verdict_from_report(report):
    if "VERDICT: ITERATE" in report:
        return "ITERATE"
    return "APPROVED"


def request_fix_code(code, report, mission, task_instruction):
    instruction = f"""
MISSION:
{mission}

TASK INSTRUCTION:
{task_instruction}

REVIEW:
{report}

ORIGINAL CODE:
{code}

Fix the code while preserving the original task objective.

Output only the corrected code.
"""
    process = subprocess.run(
        ["python", Path(__file__).parent / "Cabo.py", "-C"],
        input=instruction,
        text=True,
        capture_output=True,
    )
    return process.stdout


def request_fix_text(text, report, mission, task_instruction):
    instruction = f"""
Fix the following text based on the review.

MISSION:
{mission}

TASK INSTRUCTION:
{task_instruction}

REVIEW:
{report}

ORIGINAL TEXT:
{text}

Output only the newer version.
"""
    process = subprocess.run(
        ["python", "Cabo.py", "-T"], input=instruction, text=True, capture_output=True
    )
    return process.stdout


def build_project_context(files):
    summary = []
    for name, content in files.items():
        if name.endswith("_review.txt"):
            continue
        if name == "ALFEREZ_PLAN.txt":
            continue
        snippet = content[:800]
        summary.append(f"""
FILE: {name}
CONTENT PREVIEW:
{snippet}
""")
    return "\n".join(summary)


if __name__ == "__main__":
    process()
