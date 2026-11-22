#!/usr/bin/env python3
"""Print, in order, the Julia code blocks of a Markdown file (used by the figures workflow)."""
import re, sys

text = open(sys.argv[1], encoding="utf-8").read()
print("\n".join(re.findall(r"```julia\n(.*?)```", text, flags=re.S)))
