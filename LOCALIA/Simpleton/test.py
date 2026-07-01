def generate_enhancement(
    model="xentriom/gemma-4-12B-coder-fable5-composer2.5-v1", code="", instruction=""
):
    system = (
        "You are a proficient coder fluent in Python and Bash.\n\n"
        "If the user provides only 'BUILD' instructions, generate the complete implementation from those specifications.\n"
        "If the user provides 'ENHANCE' code, refactor it for robustness, use more idiomatic patterns, add meaningful comments at critical spots, and fix bugs while preserving intent.\\n\\n"
        "When both are provided, first enhance the existing code then implement any additions specified in the instruction block."
    )

    mode = "ENHANCE"
    if not code:
        mode = "BUILD"
    elif instruction:
        mode = "BOTH"

    prompt_segments = [system]
    if mode == "ENHANCE":
        prompt_segments.append("\n\nMODE: ENHANCE\n")
        prompt_segments.append(f"\nCODE:\n{code}")
    elif mode == "BUILD":
        prompt_segments.append("\n\nMODE: BUILD\n")
        prompt_segments.append(f"\nINSTRUCTION:\n{instruction}")
    else:
        prompt_segments.append("\n\nMODE: BOTH\n")
        prompt_segments.append(f"ENHANCED CODE:\n{code}")
        prompt_segments.append(f"\nADDITIONAL INSTRUCTIONS:\n{instruction}")

    full_prompt = "\n".join([segment for segment in prompt_segments if segment])
    return f"{model}\n{full_prompt}"
