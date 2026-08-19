# -*- coding: utf-8 -*-
"""PyMOL antibody surface hydrophobic residue analyzer.

Commands registered after loading this file:
    surface_hydrophobics
    ab_hydro

Example:
    run C:/path/pymol_antibody_surface_hydrophobic.py
    ab_hydro antibody
"""

from __future__ import print_function

import csv
import math
import os
import re
from collections import defaultdict

from pymol import cmd


# Theoretical maximum ASA values (Angstrom^2), Gly-X-Gly scale.
# Tien et al., PLoS ONE 2013, 8(11): e80635.
MAX_ASA_TIEN = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 223.0, "GLU": 225.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

DEFAULT_HYDROPHOBIC = "AVILMFWY"
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}


def _safe_name(value):
    """Return a PyMOL-friendly name."""
    value = re.sub(r"[^A-Za-z0-9_]", "_", str(value))
    return value.strip("_") or "hydro"


def _centroid(coordinates):
    """Return the arithmetic centroid of a coordinate list."""
    count = float(len(coordinates))
    return tuple(
        sum(coord[index] for coord in coordinates) / count
        for index in range(3)
    )


def _residue_selection(obj, segi, chain, resi):
    """Build a selection for one residue in the generated view object."""
    parts = [obj, "resi %s" % resi]

    if segi:
        parts.append("segi %s" % segi)

    if chain:
        parts.append("chain %s" % chain)
    else:
        parts.append('chain ""')

    return "(" + " and ".join(parts) + ")"


def surface_hydrophobics(
    selection="polymer.protein",
    context="",
    rsa_cutoff=0.25,
    sidechain_sasa_cutoff=15.0,
    hydrophobic=DEFAULT_HYDROPHOBIC,
    patch_distance=7.5,
    dot_density=3,
    state=1,
    prefix="hydro",
    csv_file="",
    include_hetatm=0,
    quiet=0,
):
    """
DESCRIPTION

    Calculate residue SASA/RSA and identify exposed hydrophobic residues.

USAGE

    surface_hydrophobics selection [, context [, rsa_cutoff
        [, sidechain_sasa_cutoff [, hydrophobic [, patch_distance
        [, dot_density [, state [, prefix [, csv_file
        [, include_hetatm [, quiet ]]]]]]]]]]]]

ARGUMENTS

    selection
        Protein residues to report.

    context
        Complete structural context used for SASA calculation.
        If empty, selection is used. For an antibody-antigen complex,
        report antibody chains in selection and pass the whole complex
        as context.

    rsa_cutoff
        Minimum residue RSA for surface exposure. Default: 0.25

    sidechain_sasa_cutoff
        Minimum side-chain SASA in Angstrom^2. Default: 15.0

    hydrophobic
        One-letter residue codes considered hydrophobic.
        Default: AVILMFWY

    patch_distance
        Maximum side-chain centroid distance used to count neighboring
        candidate residues. Default: 7.5 Angstrom

    dot_density
        PyMOL surface sampling density, 1 to 4. Default: 3

    state
        Coordinate state to analyze. Default: 1

    prefix
        Prefix for output object, selections, and default CSV.

    csv_file
        Full CSV output path. If empty, <prefix>_rsa.csv is written to
        the current PyMOL working directory.

    include_hetatm
        Keep non-protein atoms as SASA-occluding context. Default: 0
        Waters and hydrogen atoms are always removed.

OUTPUT

    <prefix>_view
        RSA-colored copy. B-factor stores RSA * 100.

    <prefix>_all_exposed
        All residues with RSA >= rsa_cutoff.

    <prefix>_exposed_hydrophobic
        Exposed hydrophobic candidates.

    <prefix>_rsa.csv
        Residue-level report.
    """

    rsa_cutoff = float(rsa_cutoff)
    sidechain_sasa_cutoff = float(sidechain_sasa_cutoff)
    patch_distance = float(patch_distance)
    dot_density = int(dot_density)
    state = int(state)
    include_hetatm = int(include_hetatm)
    quiet = int(quiet)
    prefix = _safe_name(prefix)

    hydrophobic = set(
        str(hydrophobic).upper().replace(",", "").replace(" ", "")
    )

    if not context:
        context = selection

    if not 0.0 <= rsa_cutoff <= 1.5:
        raise ValueError("rsa_cutoff should normally be between 0 and 1")

    if dot_density not in (1, 2, 3, 4):
        raise ValueError("dot_density must be 1, 2, 3, or 4")

    if state < 1:
        raise ValueError("state must be >= 1")

    if cmd.count_atoms(selection) == 0:
        raise ValueError("analysis selection contains no atoms: %s" % selection)

    if cmd.count_atoms(context) == 0:
        raise ValueError("context selection contains no atoms: %s" % context)

    # Residue identifiers to report. The context can contain additional chains.
    analyze_keys = set()
    cmd.iterate(
        "(%s) and polymer.protein" % selection,
        "analyze_keys.add((segi, chain, resi, resn))",
        space={"analyze_keys": analyze_keys},
    )

    temp_obj = cmd.get_unused_name("__rsa_context")
    view_obj = prefix + "_view"
    exposed_sel = prefix + "_all_exposed"
    candidate_sel = prefix + "_exposed_hydrophobic"

    cmd.create(temp_obj, context, state, 1, zoom=0)

    try:
        cmd.remove("%s and solvent" % temp_obj)
        cmd.remove("%s and hydro" % temp_obj)

        if not include_hetatm:
            cmd.remove("%s and not polymer.protein" % temp_obj)

        cmd.flag("ignore", temp_obj, "clear")
        cmd.set("dot_solvent", 1, temp_obj)
        cmd.set("dot_density", dot_density, temp_obj)
        cmd.get_area(temp_obj, state=1, load_b=1)

        atoms = []
        cmd.iterate_state(
            1,
            "%s and polymer.protein" % temp_obj,
            "atoms.append((segi,chain,resi,resv,resn,name,float(b),"
            "(float(x),float(y),float(z))))",
            space={"atoms": atoms},
        )

        grouped = defaultdict(list)

        for atom in atoms:
            segi, chain, resi, resv, resn, name, area, xyz = atom
            key = (segi, chain, resi, resn)

            if key not in analyze_keys:
                continue

            grouped[key].append(
                {
                    "name": name,
                    "area": area,
                    "xyz": xyz,
                    "resv": resv,
                }
            )

        rows = []

        for key, residue_atoms in grouped.items():
            segi, chain, resi, resn = key

            if resn not in MAX_ASA_TIEN:
                continue

            residue_sasa = sum(atom["area"] for atom in residue_atoms)
            sidechain_atoms = [
                atom
                for atom in residue_atoms
                if atom["name"] not in BACKBONE_NAMES
            ]
            sidechain_sasa = sum(atom["area"] for atom in sidechain_atoms)

            centroid_atoms = sidechain_atoms or residue_atoms
            centroid = _centroid(
                [atom["xyz"] for atom in centroid_atoms]
            )

            rsa = residue_sasa / MAX_ASA_TIEN[resn]
            aa = ONE_LETTER[resn]
            exposed = rsa >= rsa_cutoff
            candidate = (
                exposed
                and aa in hydrophobic
                and sidechain_sasa >= sidechain_sasa_cutoff
            )

            rows.append(
                {
                    "segi": segi,
                    "chain": chain,
                    "resi": resi,
                    "resv": residue_atoms[0]["resv"],
                    "resn": resn,
                    "aa": aa,
                    "sasa": residue_sasa,
                    "rsa": rsa,
                    "sidechain_sasa": sidechain_sasa,
                    "exposed": exposed,
                    "candidate": candidate,
                    "centroid": centroid,
                    "neighbors": 0,
                }
            )

        candidates = [row for row in rows if row["candidate"]]

        # Count candidate residues whose side-chain centroids are close.
        for index_i, row_i in enumerate(candidates):
            xi, yi, zi = row_i["centroid"]
            neighbor_count = 0

            for index_j, row_j in enumerate(candidates):
                if index_i == index_j:
                    continue

                xj, yj, zj = row_j["centroid"]
                distance = math.sqrt(
                    (xi - xj) ** 2
                    + (yi - yj) ** 2
                    + (zi - zj) ** 2
                )

                if distance <= patch_distance:
                    neighbor_count += 1

            row_i["neighbors"] = neighbor_count

        rows.sort(
            key=lambda row: (
                row["segi"],
                row["chain"],
                row["resv"],
                row["resi"],
            )
        )

        candidates.sort(
            key=lambda row: (
                -row["neighbors"],
                -row["sidechain_sasa"],
                -row["rsa"],
            )
        )

        # Keep the calculated context as the view object.
        cmd.delete(view_obj)
        cmd.set_name(temp_obj, view_obj)
        temp_obj = None

        rsa_map = {
            (row["segi"], row["chain"], row["resi"], row["resn"]):
                min(row["rsa"], 1.0) * 100.0
            for row in rows
        }

        cmd.alter(
            view_obj,
            "b=rsa_map.get((segi,chain,resi,resn),0.0)",
            space={"rsa_map": rsa_map},
        )

        exposed_expression = " or ".join(
            _residue_selection(
                view_obj, row["segi"], row["chain"], row["resi"]
            )
            for row in rows
            if row["exposed"]
        )

        candidate_expression = " or ".join(
            _residue_selection(
                view_obj, row["segi"], row["chain"], row["resi"]
            )
            for row in candidates
        )

        cmd.select(
            exposed_sel,
            exposed_expression if exposed_expression else "none",
        )
        cmd.select(
            candidate_sel,
            candidate_expression if candidate_expression else "none",
        )

        # Visualization. Blue-white-red encodes RSA, not hydrophobicity.
        cmd.hide("everything", view_obj)
        cmd.show("cartoon", view_obj)
        cmd.show("surface", view_obj)
        cmd.set("transparency", 0.35, view_obj)
        cmd.spectrum(
            "b",
            "blue_white_red",
            view_obj,
            minimum=0,
            maximum=100,
        )
        cmd.show("sticks", candidate_sel)
        cmd.color("orange", candidate_sel)

        # Use constant label expressions for compatibility with older PyMOL.
        # Do not pass the unsupported cmd.label(..., space=...) argument.
        for row in candidates:
            residue_sel = _residue_selection(
                view_obj, row["segi"], row["chain"], row["resi"]
            )
            cmd.label(
                "%s and name CA" % residue_sel,
                '"%s%s"' % (row["aa"], row["resi"]),
            )

        if not csv_file:
            csv_file = os.path.abspath(prefix + "_rsa.csv")
        else:
            csv_file = os.path.abspath(os.path.expanduser(csv_file))

        output_dir = os.path.dirname(csv_file)
        if output_dir and not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        neighbor_column = (
            "hydrophobic_neighbors_within_%.1fA" % patch_distance
        )

        fields = [
            "chain",
            "segi",
            "resi",
            "resn",
            "aa",
            "sasa_A2",
            "rsa",
            "sidechain_sasa_A2",
            "exposed",
            "hydrophobic_candidate",
            neighbor_column,
        ]

        with open(csv_file, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()

            for row in rows:
                writer.writerow(
                    {
                        "chain": row["chain"],
                        "segi": row["segi"],
                        "resi": row["resi"],
                        "resn": row["resn"],
                        "aa": row["aa"],
                        "sasa_A2": "%.3f" % row["sasa"],
                        "rsa": "%.4f" % row["rsa"],
                        "sidechain_sasa_A2":
                            "%.3f" % row["sidechain_sasa"],
                        "exposed": int(row["exposed"]),
                        "hydrophobic_candidate": int(row["candidate"]),
                        neighbor_column: row["neighbors"],
                    }
                )

        print("")
        print("Surface hydrophobic analysis finished")
        print("  Protein residues analyzed : %d" % len(rows))
        print(
            "  Exposed residues          : %d (RSA >= %.2f)"
            % (sum(row["exposed"] for row in rows), rsa_cutoff)
        )
        print("  Hydrophobic candidates    : %d" % len(candidates))
        print("  Candidate selection       : %s" % candidate_sel)
        print("  RSA-colored object        : %s" % view_obj)
        print("  CSV report                : %s" % csv_file)

        if candidates and not quiet:
            print("")
            print(
                "Ranked candidates "
                "(neighbors, side-chain SASA, RSA):"
            )

            for row in candidates:
                label = "%s/%s%s" % (
                    row["chain"] or "-",
                    row["aa"],
                    row["resi"],
                )
                print(
                    "  %-12s neighbors=%2d  scSASA=%7.2f A^2  RSA=%.3f"
                    % (
                        label,
                        row["neighbors"],
                        row["sidechain_sasa"],
                        row["rsa"],
                    )
                )

        return candidates

    finally:
        if temp_obj:
            cmd.delete(temp_obj)


def ab_hydro(
    obj="",
    chains="H+L",
    rsa_cutoff=0.25,
    sidechain_sasa_cutoff=15.0,
    prefix="ab_hydro",
    csv_file="",
    include_hetatm=0,
    quiet=0,
):
    """
DESCRIPTION

    Convenience command for antibody analysis.

USAGE

    ab_hydro [ object [, chains [, rsa_cutoff
        [, sidechain_sasa_cutoff [, prefix [, csv_file
        [, include_hetatm [, quiet ]]]]]]]]

EXAMPLES

    ab_hydro
    ab_hydro antibody
    ab_hydro antibody, chains=H+L
    ab_hydro nanobody, chains=A
    """

    if not obj:
        objects = cmd.get_object_list("enabled")

        if len(objects) == 0:
            raise ValueError("no enabled molecular object was found")

        if len(objects) > 1:
            raise ValueError(
                "multiple objects are enabled; specify one, "
                "for example: ab_hydro antibody"
            )

        obj = objects[0]

    antibody_selection = (
        "(%s) and polymer.protein and chain %s" % (obj, chains)
    )

    if cmd.count_atoms(antibody_selection) == 0:
        raise ValueError(
            "no atoms were found for chain(s) %s in %s; "
            "check the chain names" % (chains, obj)
        )

    return surface_hydrophobics(
        selection=antibody_selection,
        context=obj,
        rsa_cutoff=rsa_cutoff,
        sidechain_sasa_cutoff=sidechain_sasa_cutoff,
        prefix=prefix,
        csv_file=csv_file,
        include_hetatm=include_hetatm,
        quiet=quiet,
    )


cmd.extend("surface_hydrophobics", surface_hydrophobics)
cmd.extend("ab_hydro", ab_hydro)

print("Loaded commands: surface_hydrophobics, ab_hydro")
