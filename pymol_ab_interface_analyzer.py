# -*- coding: utf-8 -*-
"""
Auditable antigen-antibody interface interaction analyzer for PyMOL.

The command reports geometry-supported interaction candidates. It does not
calculate interaction energies and must not be interpreted as experimental
proof of a force. Rules and thresholds are written to every JSON report.

Typical use:
    run C:/path/pymol_ab_interface_analyzer.py
    ab_interface (complex and chain H+L), (complex and chain A), prefix=demo

Version: 0.1.0
License: MIT
"""

from __future__ import print_function

import csv
import json
import math
import os
from collections import defaultdict

try:
    from pymol import cmd
except ImportError:  # Allows geometry unit tests outside PyMOL.
    cmd = None


__version__ = "0.1.0"

AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

HYDROPHOBIC_RESIDUES = set(("ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR"))
BACKBONE_NAMES = set(("N", "CA", "C", "O", "OXT"))

ELEMENT_VDW = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
}

AROMATIC_RINGS = {
    "PHE": {
        "atoms": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
        "normal": ("CG", "CD1", "CD2"),
    },
    "TYR": {
        "atoms": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
        "normal": ("CG", "CD1", "CD2"),
    },
    "HIS": {
        "atoms": ("CG", "ND1", "CD2", "CE1", "NE2"),
        "normal": ("CG", "ND1", "CD2"),
    },
    # The six-membered benzene portion of indole gives a stable plane.
    "TRP": {
        "atoms": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
        "normal": ("CD2", "CE2", "CE3"),
    },
}

DEFAULTS = {
    "contact_cutoff": 4.0,
    "bsa_residue_cutoff": 1.0,
    "hbond_cutoff": 3.5,
    "hbond_angle": 63.0,
    "salt_cutoff": 4.0,
    "hydrophobic_cutoff": 4.5,
    "pi_stack_cutoff": 7.0,
    "t_stack_cutoff": 5.0,
    "pi_normal_angle": 30.0,
    "t_normal_tolerance": 30.0,
    "pi_psi_angle": 45.0,
    "cation_pi_cutoff": 6.0,
    "cation_pi_angle": 60.0,
    "vdw_tolerance": 0.5,
    "clash_overlap": 0.6,
    "disulfide_min": 1.8,
    "disulfide_max": 2.3,
}

COLORS = {
    "hydrogen_bond": "cyan",
    "salt_bridge": "magenta",
    "hydrophobic_contact": "orange",
    "pi_stacking": "purple",
    "t_stacking": "violet",
    "cation_pi": "marine",
    "disulfide": "yellow",
    "vdw_contact": "gray70",
    "steric_clash": "red",
    "close_contact": "gray50",
}


def _as_float(value, name):
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError("%s must be a number, got %r" % (name, value))


def _as_int(value, name):
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError("%s must be an integer, got %r" % (name, value))


def _truthy(value):
    return str(value).strip().lower() not in ("0", "false", "no", "off", "")


def _distance_xyz(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    n = _norm(a)
    if n < 1.0e-8:
        return None
    return (a[0] / n, a[1] / n, a[2] / n)


def _acute_angle_degrees(a, b):
    ua = _unit(a)
    ub = _unit(b)
    if ua is None or ub is None:
        return None
    cosine = max(-1.0, min(1.0, abs(_dot(ua, ub))))
    return math.degrees(math.acos(cosine))


def _centroid(coords):
    count = float(len(coords))
    return (
        sum(x[0] for x in coords) / count,
        sum(x[1] for x in coords) / count,
        sum(x[2] for x in coords) / count,
    )


def _safe_text(value):
    if value is None:
        return ""
    return str(value)


def _residue_id(atom):
    return (
        _safe_text(atom.get("segi")),
        _safe_text(atom.get("chain")),
        _safe_text(atom.get("resi")),
        _safe_text(atom.get("resn")).upper(),
    )


def _residue_key_text(residue_id):
    segi, chain, resi, resn = residue_id
    prefix = ""
    if segi:
        prefix += segi + "/"
    prefix += (chain if chain else "_")
    return "%s/%s%s" % (prefix, resn, resi)


def _residue_label(residue_id):
    resn = residue_id[3]
    return "%s%s" % (AA1.get(resn, "X"), residue_id[2])


def _atom_label(atom):
    return "%s:%s" % (_residue_key_text(_residue_id(atom)), atom["name"])


def _atom_match_key(atom):
    return (
        _safe_text(atom.get("segi")),
        _safe_text(atom.get("chain")),
        _safe_text(atom.get("resi")),
        _safe_text(atom.get("resn")).upper(),
        _safe_text(atom.get("name")).upper(),
        _safe_text(atom.get("alt")),
    )


def _atom_uid(atom):
    return (_safe_text(atom.get("model")), int(atom["index"]))


def _element_vdw(atom):
    value = atom.get("vdw", 0.0)
    try:
        value = float(value)
    except (TypeError, ValueError):
        value = 0.0
    if value > 0.1:
        return value
    return ELEMENT_VDW.get(_safe_text(atom.get("elem")).upper(), 1.70)


def _get_atoms(selection, state):
    model = cmd.get_model("(%s) and not hydro and not solvent" % selection, state)
    atoms = []
    for atom in model.atom:
        atoms.append({
            "model": _safe_text(getattr(atom, "model", "")),
            "index": int(atom.index),
            "segi": _safe_text(getattr(atom, "segi", "")),
            "chain": _safe_text(getattr(atom, "chain", "")),
            "resi": _safe_text(getattr(atom, "resi", "")),
            "resn": _safe_text(getattr(atom, "resn", "")).upper(),
            "name": _safe_text(getattr(atom, "name", "")).upper(),
            "alt": _safe_text(getattr(atom, "alt", "")),
            "elem": _safe_text(getattr(atom, "symbol", getattr(atom, "elem", ""))).upper(),
            "coord": tuple(float(x) for x in atom.coord),
            "vdw": float(getattr(atom, "vdw", 0.0) or 0.0),
            "b": float(getattr(atom, "b", 0.0) or 0.0),
        })
    return atoms


def _grid_pairs(atoms1, atoms2, cutoff):
    """Return all atom pairs within cutoff using a dependency-free cell list."""
    cutoff = float(cutoff)
    if cutoff <= 0.0:
        return []
    grid = defaultdict(list)
    inv = 1.0 / cutoff
    for atom in atoms2:
        x, y, z = atom["coord"]
        cell = (int(math.floor(x * inv)), int(math.floor(y * inv)), int(math.floor(z * inv)))
        grid[cell].append(atom)

    cutoff2 = cutoff * cutoff
    pairs = []
    for atom1 in atoms1:
        x, y, z = atom1["coord"]
        base = (int(math.floor(x * inv)), int(math.floor(y * inv)), int(math.floor(z * inv)))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for atom2 in grid.get((base[0] + dx, base[1] + dy, base[2] + dz), ()):
                        ax = x - atom2["coord"][0]
                        ay = y - atom2["coord"][1]
                        az = z - atom2["coord"][2]
                        d2 = ax * ax + ay * ay + az * az
                        if d2 <= cutoff2:
                            pairs.append((atom1, atom2, math.sqrt(d2)))
    return pairs


def _best_by_residue_pair(atom_pairs):
    best = {}
    counts = defaultdict(int)
    for atom1, atom2, distance in atom_pairs:
        key = (_residue_id(atom1), _residue_id(atom2))
        counts[key] += 1
        if key not in best or distance < best[key][2]:
            best[key] = (atom1, atom2, distance)
    return best, counts


def _add_interaction(store, interaction_type, atom1, atom2, distance,
                     confidence, criterion, geometry=""):
    r1 = _residue_id(atom1)
    r2 = _residue_id(atom2)
    key = (interaction_type, r1, r2)
    row = {
        "type": interaction_type,
        "confidence": confidence,
        "partner1_residue": _residue_key_text(r1),
        "partner1_atom": atom1["name"],
        "partner2_residue": _residue_key_text(r2),
        "partner2_atom": atom2["name"],
        "distance_A": round(float(distance), 3),
        "geometry": geometry,
        "criterion": criterion,
        "_atom1": atom1,
        "_atom2": atom2,
        "_r1": r1,
        "_r2": r2,
    }
    if key not in store or distance < store[key]["distance_A"]:
        store[key] = row


def _ring_records(atoms, side, warnings):
    grouped = defaultdict(dict)
    representative = {}
    for atom in atoms:
        rid = _residue_id(atom)
        grouped[rid][atom["name"]] = atom
        representative[rid] = atom

    rings = []
    for rid, names in grouped.items():
        definition = AROMATIC_RINGS.get(rid[3])
        if definition is None:
            continue
        missing = [name for name in definition["atoms"] if name not in names]
        if missing:
            warnings.append(
                "%s aromatic residue %s lacks ring atoms: %s"
                % (side, _residue_key_text(rid), ",".join(missing))
            )
            continue
        coords = [names[name]["coord"] for name in definition["atoms"]]
        n0, n1, n2 = (names[name]["coord"] for name in definition["normal"])
        normal = _unit(_cross(_sub(n1, n0), _sub(n2, n0)))
        if normal is None:
            warnings.append("%s ring %s is geometrically degenerate" % (side, _residue_key_text(rid)))
            continue
        rings.append({
            "residue": rid,
            "centroid": _centroid(coords),
            "normal": normal,
            "atom": names[definition["normal"][0]],
        })
    return rings


def _cation_records(atoms):
    grouped = defaultdict(dict)
    for atom in atoms:
        grouped[_residue_id(atom)][atom["name"]] = atom

    records = []
    for rid, names in grouped.items():
        if rid[3] == "LYS" and "NZ" in names:
            records.append({"residue": rid, "center": names["NZ"]["coord"], "atom": names["NZ"]})
        elif rid[3] == "ARG" and all(name in names for name in ("CZ", "NH1", "NH2")):
            records.append({
                "residue": rid,
                "center": _centroid([names["CZ"]["coord"], names["NH1"]["coord"], names["NH2"]["coord"]]),
                "atom": names["CZ"],
            })
    return records


def _ring_interactions(rings1, rings2, store, cfg):
    for ring1 in rings1:
        for ring2 in rings2:
            vector = _sub(ring2["centroid"], ring1["centroid"])
            distance = _norm(vector)
            normal_angle = _acute_angle_degrees(ring1["normal"], ring2["normal"])
            psi1 = _acute_angle_degrees(ring1["normal"], vector)
            psi2 = _acute_angle_degrees(ring2["normal"], vector)
            if normal_angle is None or psi1 is None or psi2 is None:
                continue
            psi = min(psi1, psi2)
            geometry = "normal_angle=%.1f deg; psi_min=%.1f deg" % (normal_angle, psi)
            if (
                distance <= cfg["pi_stack_cutoff"]
                and normal_angle <= cfg["pi_normal_angle"]
                and psi <= cfg["pi_psi_angle"]
            ):
                _add_interaction(
                    store, "pi_stacking", ring1["atom"], ring2["atom"], distance,
                    "geometry_supported",
                    "centroid<=%.1f A; normal_angle<=%.1f deg; psi<=%.1f deg"
                    % (cfg["pi_stack_cutoff"], cfg["pi_normal_angle"], cfg["pi_psi_angle"]),
                    geometry,
                )
            elif (
                distance <= cfg["t_stack_cutoff"]
                and abs(normal_angle - 90.0) <= cfg["t_normal_tolerance"]
                and psi <= cfg["pi_psi_angle"]
            ):
                _add_interaction(
                    store, "t_stacking", ring1["atom"], ring2["atom"], distance,
                    "geometry_supported",
                    "centroid<=%.1f A; |normal_angle-90|<=%.1f deg; psi<=%.1f deg"
                    % (cfg["t_stack_cutoff"], cfg["t_normal_tolerance"], cfg["pi_psi_angle"]),
                    geometry,
                )


def _cation_pi_interactions(cations, rings, cation_is_partner1, store, cfg):
    for cation in cations:
        for ring in rings:
            vector = _sub(cation["center"], ring["centroid"])
            distance = _norm(vector)
            angle = _acute_angle_degrees(ring["normal"], vector)
            if angle is None:
                continue
            if distance <= cfg["cation_pi_cutoff"] and angle <= cfg["cation_pi_angle"]:
                geometry = "ring_axis_angle=%.1f deg" % angle
                criterion = "cation-ring centroid<=%.1f A; ring_axis_angle<=%.1f deg" % (
                    cfg["cation_pi_cutoff"], cfg["cation_pi_angle"]
                )
                if cation_is_partner1:
                    _add_interaction(
                        store, "cation_pi", cation["atom"], ring["atom"], distance,
                        "geometry_supported; charge_state_assumed", criterion, geometry,
                    )
                else:
                    _add_interaction(
                        store, "cation_pi", ring["atom"], cation["atom"], distance,
                        "geometry_supported; charge_state_assumed", criterion, geometry,
                    )


def _py_quote(value):
    return '"%s"' % _safe_text(value).replace("\\", "\\\\").replace('"', '\\"')


def _residue_expression(base_selection, residue_ids):
    parts = []
    for segi, chain, resi, resn in sorted(residue_ids):
        clauses = ["(%s)" % base_selection, "resi %s" % _py_quote(resi), "resn %s" % resn]
        clauses.append("chain %s" % _py_quote(chain))
        if segi:
            clauses.append("segi %s" % _py_quote(segi))
        parts.append("(" + " and ".join(clauses) + ")")
    return "none" if not parts else " or ".join(parts)


def _atom_expression(atom):
    model = atom["model"]
    if model:
        return "(model %s and index %d)" % (_py_quote(model), atom["index"])
    return "(index %d)" % atom["index"]


def _validate_unique_match_keys(atoms1, atoms2):
    keys1 = [_atom_match_key(atom) for atom in atoms1]
    keys2 = [_atom_match_key(atom) for atom in atoms2]
    duplicate1 = len(keys1) != len(set(keys1))
    duplicate2 = len(keys2) != len(set(keys2))
    overlap = set(keys1).intersection(keys2)
    if duplicate1 or duplicate2 or overlap:
        raise ValueError(
            "Atom identifiers are not unique across the two partners. "
            "Use one complex object with distinct chain IDs (and distinct segi IDs if needed)."
        )


def _calculate_bsa(selection1, selection2, state, residue_cutoff, warnings):
    temp1 = cmd.get_unused_name("_abint_p1")
    temp2 = cmd.get_unused_name("_abint_p2")
    temp_complex = cmd.get_unused_name("_abint_complex")
    settings = {}
    try:
        for setting in ("dot_solvent", "dot_density"):
            try:
                settings[setting] = cmd.get(setting)
            except Exception:
                pass

        cmd.create(temp1, "(%s) and not hydro and not solvent" % selection1, state, 1)
        cmd.create(temp2, "(%s) and not hydro and not solvent" % selection2, state, 1)
        cmd.create(temp_complex, "%s or %s" % (temp1, temp2), 1, 1)

        cmd.set("dot_solvent", 1)
        cmd.set("dot_density", 4)

        cmd.get_area(temp1, state=1, load_b=1)
        isolated1 = {_atom_match_key(a): a["b"] for a in _get_atoms(temp1, 1)}
        cmd.get_area(temp2, state=1, load_b=1)
        isolated2 = {_atom_match_key(a): a["b"] for a in _get_atoms(temp2, 1)}
        cmd.get_area(temp_complex, state=1, load_b=1)
        complex_area = {_atom_match_key(a): a["b"] for a in _get_atoms(temp_complex, 1)}

        delta_by_residue = defaultdict(float)
        total1 = 0.0
        total2 = 0.0

        for atom in _get_atoms(temp1, 1):
            key = _atom_match_key(atom)
            if key not in complex_area:
                warnings.append("SASA key missing in complex: %s" % (key,))
                continue
            delta = max(0.0, isolated1.get(key, 0.0) - complex_area[key])
            delta_by_residue[("partner1", _residue_id(atom))] += delta
            total1 += delta

        for atom in _get_atoms(temp2, 1):
            key = _atom_match_key(atom)
            if key not in complex_area:
                warnings.append("SASA key missing in complex: %s" % (key,))
                continue
            delta = max(0.0, isolated2.get(key, 0.0) - complex_area[key])
            delta_by_residue[("partner2", _residue_id(atom))] += delta
            total2 += delta

        interface = {
            side: set(rid for (which, rid), delta in delta_by_residue.items()
                      if which == side and delta >= residue_cutoff)
            for side in ("partner1", "partner2")
        }
        return {
            "delta_by_residue": delta_by_residue,
            "interface": interface,
            "buried_area_partner1_A2": total1,
            "buried_area_partner2_A2": total2,
            "bsa_A2": 0.5 * (total1 + total2),
            "side_asymmetry_A2": abs(total1 - total2),
        }
    finally:
        for name in (temp1, temp2, temp_complex):
            try:
                cmd.delete(name)
            except Exception:
                pass
        for setting, value in settings.items():
            try:
                cmd.set(setting, value)
            except Exception:
                pass


def _find_hbonds(selection1, selection2, state, atoms1, atoms2, store, cfg, warnings):
    lookup1 = {_atom_uid(atom): atom for atom in atoms1}
    lookup2 = {_atom_uid(atom): atom for atom in atoms2}
    calls = (
        (
            "(%s) and donors" % selection1,
            "(%s) and acceptors" % selection2,
            lookup1,
            lookup2,
            False,
        ),
        (
            "(%s) and acceptors" % selection1,
            "(%s) and donors" % selection2,
            lookup1,
            lookup2,
            False,
        ),
    )
    found = 0
    for sel1, sel2, left_lookup, right_lookup, reverse in calls:
        try:
            pairs = cmd.find_pairs(
                sel1, sel2, state1=state, state2=state,
                cutoff=cfg["hbond_cutoff"], mode=1, angle=cfg["hbond_angle"]
            )
        except Exception as exc:
            warnings.append("PyMOL hydrogen-bond search failed: %s" % exc)
            return
        for uid1, uid2 in pairs:
            atom1 = left_lookup.get((_safe_text(uid1[0]), int(uid1[1])))
            atom2 = right_lookup.get((_safe_text(uid2[0]), int(uid2[1])))
            if atom1 is None or atom2 is None:
                continue
            distance = _distance_xyz(atom1["coord"], atom2["coord"])
            _add_interaction(
                store, "hydrogen_bond", atom1, atom2, distance,
                "putative; PyMOL donor/acceptor and orientation",
                "donor-acceptor<=%.1f A; PyMOL mode=1; angle<=%.1f deg"
                % (cfg["hbond_cutoff"], cfg["hbond_angle"]),
            )
            found += 1
    if found == 0:
        warnings.append(
            "No putative hydrogen bonds were found. Check protonation, atom typing, "
            "missing side-chain atoms, and structure resolution."
        )


def _interaction_rows(atoms1, atoms2, selection1, selection2, state, cfg, warnings):
    store = {}

    max_cutoff = max(
        cfg["contact_cutoff"], cfg["hydrophobic_cutoff"],
        cfg["salt_cutoff"], cfg["disulfide_max"],
        max(ELEMENT_VDW.values()) * 2.0 + cfg["vdw_tolerance"],
    )
    nearby = _grid_pairs(atoms1, atoms2, max_cutoff)

    contact_pairs = [(a, b, d) for a, b, d in nearby if d <= cfg["contact_cutoff"]]
    best_contacts, contact_counts = _best_by_residue_pair(contact_pairs)
    for key, (atom1, atom2, distance) in best_contacts.items():
        _add_interaction(
            store, "close_contact", atom1, atom2, distance,
            "distance_defined",
            "any heavy-atom distance<=%.1f A" % cfg["contact_cutoff"],
            "atom_pair_count=%d" % contact_counts[key],
        )

    positive = set((("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")))
    negative = set((("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")))
    for atom1, atom2, distance in nearby:
        type1 = (atom1["resn"], atom1["name"])
        type2 = (atom2["resn"], atom2["name"])
        if distance <= cfg["salt_cutoff"] and (
            (type1 in positive and type2 in negative)
            or (type1 in negative and type2 in positive)
        ):
            _add_interaction(
                store, "salt_bridge", atom1, atom2, distance,
                "geometry_supported; standard side-chain charge assumed",
                "oppositely charged side-chain atoms<=%.1f A; HIS/termini excluded"
                % cfg["salt_cutoff"],
            )

        hydrophobic1 = (
            atom1["resn"] in HYDROPHOBIC_RESIDUES
            and atom1["name"] not in BACKBONE_NAMES
            and atom1["elem"] in ("C", "S")
        )
        hydrophobic2 = (
            atom2["resn"] in HYDROPHOBIC_RESIDUES
            and atom2["name"] not in BACKBONE_NAMES
            and atom2["elem"] in ("C", "S")
        )
        if hydrophobic1 and hydrophobic2 and distance <= cfg["hydrophobic_cutoff"]:
            _add_interaction(
                store, "hydrophobic_contact", atom1, atom2, distance,
                "distance_defined",
                "hydrophobic side-chain C/S atoms<=%.1f A" % cfg["hydrophobic_cutoff"],
            )

        if (
            atom1["resn"] == "CYS" and atom1["name"] == "SG"
            and atom2["resn"] == "CYS" and atom2["name"] == "SG"
            and cfg["disulfide_min"] <= distance <= cfg["disulfide_max"]
        ):
            _add_interaction(
                store, "disulfide", atom1, atom2, distance,
                "strong_geometry",
                "SG-SG %.1f-%.1f A; bond order not assigned"
                % (cfg["disulfide_min"], cfg["disulfide_max"]),
            )

        vdw_gap = distance - (_element_vdw(atom1) + _element_vdw(atom2))
        if vdw_gap < -cfg["clash_overlap"]:
            _add_interaction(
                store, "steric_clash", atom1, atom2, distance,
                "geometry_warning",
                "distance < vdW-radius sum - %.1f A" % cfg["clash_overlap"],
                "vdw_gap=%.3f A" % vdw_gap,
            )
        elif abs(vdw_gap) <= cfg["vdw_tolerance"]:
            _add_interaction(
                store, "vdw_contact", atom1, atom2, distance,
                "radius_defined",
                "|distance-vdW-radius sum|<=%.1f A" % cfg["vdw_tolerance"],
                "vdw_gap=%.3f A" % vdw_gap,
            )

    _find_hbonds(selection1, selection2, state, atoms1, atoms2, store, cfg, warnings)

    rings1 = _ring_records(atoms1, "partner1", warnings)
    rings2 = _ring_records(atoms2, "partner2", warnings)
    _ring_interactions(rings1, rings2, store, cfg)

    cations1 = _cation_records(atoms1)
    cations2 = _cation_records(atoms2)
    _cation_pi_interactions(cations1, rings2, True, store, cfg)
    _cation_pi_interactions(cations2, rings1, False, store, cfg)

    rows = list(store.values())
    rows.sort(key=lambda row: (
        row["type"], row["partner1_residue"], row["partner2_residue"], row["distance_A"]
    ))
    return rows, best_contacts, contact_counts


def _make_visuals(prefix, selection1, selection2, rows, interface1, interface2, labels):
    names = {}
    names["p1_interface"] = prefix + "_p1_interface"
    names["p2_interface"] = prefix + "_p2_interface"
    names["all_interface"] = prefix + "_interface"

    cmd.select(names["p1_interface"], _residue_expression(selection1, interface1))
    cmd.select(names["p2_interface"], _residue_expression(selection2, interface2))
    cmd.select(
        names["all_interface"],
        "%s or %s" % (names["p1_interface"], names["p2_interface"]),
    )

    cmd.hide("everything", "(%s) or (%s)" % (selection1, selection2))
    cmd.show("cartoon", "(%s) or (%s)" % (selection1, selection2))
    cmd.color("lightblue", selection1)
    cmd.color("lightpink", selection2)
    cmd.show("sticks", names["all_interface"])
    cmd.color("cyan", names["p1_interface"])
    cmd.color("salmon", names["p2_interface"])

    distance_names = []
    for interaction_type in COLORS:
        typed_rows = [row for row in rows if row["type"] == interaction_type]
        if not typed_rows:
            continue
        obj = prefix + "_" + interaction_type
        for row in typed_rows:
            try:
                cmd.distance(obj, _atom_expression(row["_atom1"]), _atom_expression(row["_atom2"]))
            except Exception:
                continue
        try:
            cmd.color(COLORS[interaction_type], obj)
            cmd.hide("labels", obj)
            cmd.set("dash_width", 2.5, obj)
        except Exception:
            pass
        distance_names.append(obj)

    if labels:
        for side_name in (names["p1_interface"], names["p2_interface"]):
            try:
                cmd.label(
                    "(%s) and name CA" % side_name,
                    '"%s%s" % (one_letter[resn], resi)',
                    space={"one_letter": AA1},
                )
            except TypeError:
                # Compatible fallback for older PyMOL versions without label(space=...).
                for rid in (interface1 if side_name == names["p1_interface"] else interface2):
                    expr = _residue_expression(side_name, set((rid,)))
                    cmd.label("(%s) and name CA" % expr, _py_quote(_residue_label(rid)))
            except Exception:
                pass

    group_members = list(names.values()) + distance_names
    try:
        cmd.group(prefix, " ".join(group_members))
    except Exception:
        pass
    try:
        cmd.orient(names["all_interface"])
        cmd.zoom(names["all_interface"], 8)
    except Exception:
        pass
    return names


def _clean_rows(rows):
    cleaned = []
    for row in rows:
        cleaned.append({key: value for key, value in row.items() if not key.startswith("_")})
    return cleaned


def _write_outputs(prefix, output_dir, rows, residue_rows, summary):
    output_dir = os.path.abspath(os.path.expanduser(output_dir or os.getcwd()))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    interactions_path = os.path.join(output_dir, prefix + "_interactions.csv")
    residues_path = os.path.join(output_dir, prefix + "_residues.csv")
    summary_path = os.path.join(output_dir, prefix + "_summary.json")

    interaction_fields = (
        "type", "confidence", "partner1_residue", "partner1_atom",
        "partner2_residue", "partner2_atom", "distance_A", "geometry", "criterion",
    )
    with open(interactions_path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=interaction_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(_clean_rows(rows))

    residue_fields = (
        "partner", "residue", "one_letter", "contact_interface",
        "bsa_interface", "delta_sasa_A2", "minimum_contact_A",
        "contact_atom_pair_count", "interaction_types",
    )
    with open(residues_path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=residue_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(residue_rows)

    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)

    return interactions_path, residues_path, summary_path


def _residue_table(atoms1, atoms2, rows, best_contacts, contact_counts, bsa, cfg):
    all_residues = {
        "partner1": set(_residue_id(atom) for atom in atoms1),
        "partner2": set(_residue_id(atom) for atom in atoms2),
    }
    contact_min = {}
    for pair, best in best_contacts.items():
        contact_min[("partner1", pair[0])] = min(
            contact_min.get(("partner1", pair[0]), 999.0), best[2]
        )
        contact_min[("partner2", pair[1])] = min(
            contact_min.get(("partner2", pair[1]), 999.0), best[2]
        )

    contact_pair_count = defaultdict(int)
    for (r1, r2), count in contact_counts.items():
        contact_pair_count[("partner1", r1)] += count
        contact_pair_count[("partner2", r2)] += count

    interaction_types = defaultdict(set)
    for row in rows:
        if row["type"] != "close_contact":
            interaction_types[("partner1", row["_r1"])].add(row["type"])
            interaction_types[("partner2", row["_r2"])].add(row["type"])

    contact1 = set(pair[0] for pair in best_contacts)
    contact2 = set(pair[1] for pair in best_contacts)
    table = []
    for side in ("partner1", "partner2"):
        for rid in sorted(all_residues[side]):
            delta = bsa["delta_by_residue"].get((side, rid), 0.0)
            is_contact = rid in (contact1 if side == "partner1" else contact2)
            is_bsa = delta >= cfg["bsa_residue_cutoff"]
            if not is_contact and not is_bsa:
                continue
            minimum = contact_min.get((side, rid))
            table.append({
                "partner": side,
                "residue": _residue_key_text(rid),
                "one_letter": _residue_label(rid),
                "contact_interface": int(is_contact),
                "bsa_interface": int(is_bsa),
                "delta_sasa_A2": round(delta, 3),
                "minimum_contact_A": "" if minimum is None else round(minimum, 3),
                "contact_atom_pair_count": contact_pair_count.get((side, rid), 0),
                "interaction_types": ";".join(sorted(interaction_types.get((side, rid), set()))),
            })
    return table, contact1, contact2


def ab_interface(
    partner1,
    partner2,
    prefix="abint",
    state=1,
    output_dir="",
    contact_cutoff=4.0,
    bsa_cutoff=1.0,
    hbond_cutoff=3.5,
    salt_cutoff=4.0,
    hydrophobic_cutoff=4.5,
    labels=1,
    quiet=0,
):
    """
    Analyze an antigen-antibody or other protein-protein interface.

    partner1, partner2
        Two non-overlapping PyMOL selections. Conventionally partner1 is the
        antibody and partner2 is the antigen, but the calculation is symmetric.

    Example
        ab_interface (complex and chain H+L), (complex and chain A), prefix=case1
    """
    if cmd is None:
        raise RuntimeError("This command must be run inside PyMOL.")

    prefix = str(prefix).strip()
    if not prefix:
        raise ValueError("prefix cannot be empty")
    state = _as_int(state, "state")
    if state < 1:
        raise ValueError("state must be >= 1")

    cfg = dict(DEFAULTS)
    cfg["contact_cutoff"] = _as_float(contact_cutoff, "contact_cutoff")
    cfg["bsa_residue_cutoff"] = _as_float(bsa_cutoff, "bsa_cutoff")
    cfg["hbond_cutoff"] = _as_float(hbond_cutoff, "hbond_cutoff")
    cfg["salt_cutoff"] = _as_float(salt_cutoff, "salt_cutoff")
    cfg["hydrophobic_cutoff"] = _as_float(hydrophobic_cutoff, "hydrophobic_cutoff")

    count1 = cmd.count_atoms("(%s) and not hydro and not solvent" % partner1)
    count2 = cmd.count_atoms("(%s) and not hydro and not solvent" % partner2)
    overlap = cmd.count_atoms("(%s) and (%s)" % (partner1, partner2))
    if count1 == 0:
        raise ValueError("partner1 selection is empty")
    if count2 == 0:
        raise ValueError("partner2 selection is empty")
    if overlap:
        raise ValueError("partner1 and partner2 overlap by %d atoms" % overlap)

    warnings = []
    atoms1 = _get_atoms(partner1, state)
    atoms2 = _get_atoms(partner2, state)
    if not atoms1 or not atoms2:
        raise ValueError("One partner has no heavy non-solvent atoms in state %d" % state)
    _validate_unique_match_keys(atoms1, atoms2)

    altlocs = sorted(set(
        atom["alt"] for atom in atoms1 + atoms2 if atom["alt"]
    ))
    if altlocs:
        warnings.append(
            "Alternate conformations are present (%s). Results include every loaded altloc; "
            "choose one conformer before analysis for an unambiguous report."
            % ",".join(altlocs)
        )

    hydrogen_count = cmd.count_atoms("((%s) or (%s)) and hydro" % (partner1, partner2))
    if hydrogen_count == 0:
        warnings.append(
            "No explicit hydrogens are loaded. Hydrogen bonds are putative and use "
            "PyMOL donor/acceptor typing plus heavy-atom orientation."
        )

    bsa = _calculate_bsa(
        partner1, partner2, state, cfg["bsa_residue_cutoff"], warnings
    )
    rows, best_contacts, contact_counts = _interaction_rows(
        atoms1, atoms2, partner1, partner2, state, cfg, warnings
    )
    residue_rows, contact1, contact2 = _residue_table(
        atoms1, atoms2, rows, best_contacts, contact_counts, bsa, cfg
    )
    combined1 = contact1.union(bsa["interface"]["partner1"])
    combined2 = contact2.union(bsa["interface"]["partner2"])

    visual_names = _make_visuals(
        prefix, partner1, partner2, rows, combined1, combined2, _truthy(labels)
    )

    type_counts = defaultdict(int)
    for row in rows:
        type_counts[row["type"]] += 1

    summary = {
        "tool": "PyMOL Antibody-Antigen Interface Analyzer",
        "version": __version__,
        "interpretation": (
            "Geometry-supported interaction candidates from a static structure; "
            "not interaction energies and not experimental proof of forces."
        ),
        "partner1_selection": partner1,
        "partner2_selection": partner2,
        "state": state,
        "atom_counts": {"partner1": len(atoms1), "partner2": len(atoms2)},
        "interface_residue_counts": {
            "partner1_contact": len(contact1),
            "partner2_contact": len(contact2),
            "partner1_bsa": len(bsa["interface"]["partner1"]),
            "partner2_bsa": len(bsa["interface"]["partner2"]),
            "partner1_union": len(combined1),
            "partner2_union": len(combined2),
        },
        "buried_surface_area": {
            "BSA_A2": round(bsa["bsa_A2"], 3),
            "partner1_buried_A2": round(bsa["buried_area_partner1_A2"], 3),
            "partner2_buried_A2": round(bsa["buried_area_partner2_A2"], 3),
            "side_asymmetry_A2": round(bsa["side_asymmetry_A2"], 3),
            "formula": "(SASA_partner1 + SASA_partner2 - SASA_complex) / 2",
            "pymol_settings": {"dot_solvent": 1, "dot_density": 4},
        },
        "interaction_counts": dict(sorted(type_counts.items())),
        "thresholds": cfg,
        "assumptions": [
            "Standard protonation is assumed for Lys, Arg, Asp and Glu side chains.",
            "Histidine and chain termini are excluded from definite salt bridges.",
            "Hydrogen bonds depend on PyMOL atom typing and orientation.",
            "Aromatic interactions use ring centroid and plane geometry.",
            "Water-mediated contacts, metal coordination and binding free energy are not assigned.",
        ],
        "warnings": warnings,
        "visual_objects": visual_names,
        "method_references": [
            {
                "name": "BindCraft",
                "use": "4.0 A heavy-atom interface contact convention; Rosetta is not embedded here.",
                "url": "https://github.com/martinpacesa/BindCraft",
            },
            {
                "name": "PDBe Arpeggio",
                "use": "Auditable atom-contact classification and explicit data-quality caveats.",
                "url": "https://github.com/PDBeurope/arpeggio",
            },
            {
                "name": "GetContacts",
                "use": "Published default geometric cutoffs for salt, H-bond and aromatic contacts.",
                "url": "https://github.com/getcontacts/getcontacts",
            },
            {
                "name": "PRODIGY",
                "use": "Protein-protein contact and buried-surface-area validation context.",
                "url": "https://github.com/haddocking/prodigy",
            },
        ],
    }

    interactions_path, residues_path, summary_path = _write_outputs(
        prefix, output_dir, rows, residue_rows, summary
    )

    if not _truthy(quiet):
        print("")
        print("=== Antibody-Antigen Interface Analysis ===")
        print("Partner 1 heavy atoms: %d" % len(atoms1))
        print("Partner 2 heavy atoms: %d" % len(atoms2))
        print(
            "Interface residues (contact/BSA/union): P1 %d/%d/%d; P2 %d/%d/%d"
            % (
                len(contact1), len(bsa["interface"]["partner1"]), len(combined1),
                len(contact2), len(bsa["interface"]["partner2"]), len(combined2),
            )
        )
        print("Buried surface area: %.1f A^2" % bsa["bsa_A2"])
        print("Interaction residue-pair counts:")
        for interaction_type in sorted(type_counts):
            print("  %-22s %d" % (interaction_type, type_counts[interaction_type]))
        if warnings:
            print("Warnings:")
            for warning in warnings:
                print("  - %s" % warning)
        print("Interactions CSV: %s" % interactions_path)
        print("Residues CSV:     %s" % residues_path)
        print("Summary JSON:     %s" % summary_path)
        print(
            "Interpretation: geometry-supported candidates, not measured forces or energies."
        )
        print("")

    return {
        "interactions": _clean_rows(rows),
        "residues": residue_rows,
        "summary": summary,
        "files": {
            "interactions": interactions_path,
            "residues": residues_path,
            "summary": summary_path,
        },
    }


def ab_interface_rules():
    """Print the default interaction definitions used by ab_interface."""
    print(json.dumps({
        "version": __version__,
        "thresholds": DEFAULTS,
        "hydrophobic_residues": sorted(HYDROPHOBIC_RESIDUES),
        "aromatic_residues": sorted(AROMATIC_RINGS),
        "charged_side_chains": {
            "positive": ["LYS", "ARG"],
            "negative": ["ASP", "GLU"],
            "excluded_as_ambiguous": ["HIS", "N-terminus", "C-terminus"],
        },
    }, indent=2, sort_keys=True))


if cmd is not None:
    cmd.extend("ab_interface", ab_interface)
    cmd.extend("ab_interface_rules", ab_interface_rules)
