"""Headless PyMOL integration validation on the 1MLC Fab-lysozyme complex."""

import json
import os

from pymol import cmd

from pymol_ab_interface_analyzer import ab_interface


cmd.fetch("1mlc", async_=0)
result = ab_interface(
    "1mlc and chain A+B",
    "1mlc and chain E",
    prefix="validation_1mlc",
    output_dir=".",
    labels=0,
    quiet=1,
)

summary = result["summary"]
counts = summary["interface_residue_counts"]
bsa = summary["buried_surface_area"]["BSA_A2"]
types = summary["interaction_counts"]

assert counts["partner1_union"] >= 8, counts
assert counts["partner2_union"] >= 8, counts
assert 300.0 <= bsa <= 2000.0, bsa
assert types.get("close_contact", 0) >= 8, types
assert types.get("hydrogen_bond", 0) >= 1, types

for path in result["files"].values():
    assert os.path.isfile(path), path

with open(result["files"]["summary"], encoding="utf-8") as handle:
    reloaded = json.load(handle)
assert reloaded["version"] == result["summary"]["version"]

print("1MLC validation passed")
print(json.dumps({
    "BSA_A2": bsa,
    "interface_residue_counts": counts,
    "interaction_counts": types,
}, indent=2, sort_keys=True))
