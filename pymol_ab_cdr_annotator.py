# -*- coding: utf-8 -*-
"""PyMOL antibody / nanobody Kabat CDR annotator.

Creates CDR selections and colors them. Default scheme is Kabat.
Numbering prefers ANARCI when HMMER is available; otherwise uses a
built-in consensus alignment that needs no external binaries.

Commands:
    ab_cdr
    ab_cdr_help

Examples:
    run C:/path/pymol_ab_cdr_annotator.py
    ab_cdr
    ab_cdr antibody
    ab_cdr H+L
"""

from __future__ import print_function

import re
from collections import defaultdict

try:
    from pymol import cmd
except ImportError:
    cmd = None


__version__ = "0.2.3"

_KABAT_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "ASX": "B", "GLX": "Z", "UNK": "X", "MSE": "M", "SEC": "U",
}

KABAT_CDR_RANGES = {
    "H": {
        "HCDR1": (31, 35),
        "HCDR2": (50, 65),
        "HCDR3": (95, 102),
    },
    "L": {
        "LCDR1": (24, 34),
        "LCDR2": (50, 56),
        "LCDR3": (89, 97),
    },
}

CDR_COLORS = {
    "HCDR1": "yellow",
    "HCDR2": "orange",
    "HCDR3": "red",
    "LCDR1": "greencyan",
    "LCDR2": "green",
    "LCDR3": "marine",
}

_BLOSUM62 = {
    ("A", "A"): 4, ("A", "R"): -1, ("A", "N"): -2, ("A", "D"): -2, ("A", "C"): 0,
    ("A", "Q"): -1, ("A", "E"): -1, ("A", "G"): 0, ("A", "H"): -2, ("A", "I"): -1,
    ("A", "L"): -1, ("A", "K"): -1, ("A", "M"): -1, ("A", "F"): -2, ("A", "P"): -1,
    ("A", "S"): 1, ("A", "T"): 0, ("A", "W"): -3, ("A", "Y"): -2, ("A", "V"): 0,
    ("R", "R"): 5, ("R", "N"): 0, ("R", "D"): -2, ("R", "C"): -3, ("R", "Q"): 1,
    ("R", "E"): 0, ("R", "G"): -2, ("R", "H"): 0, ("R", "I"): -3, ("R", "L"): -2,
    ("R", "K"): 2, ("R", "M"): -1, ("R", "F"): -3, ("R", "P"): -2, ("R", "S"): -1,
    ("R", "T"): -1, ("R", "W"): -3, ("R", "Y"): -2, ("R", "V"): -3,
    ("N", "N"): 6, ("N", "D"): 1, ("N", "C"): -3, ("N", "Q"): 0, ("N", "E"): 0,
    ("N", "G"): 0, ("N", "H"): 1, ("N", "I"): -3, ("N", "L"): -3, ("N", "K"): 0,
    ("N", "M"): -2, ("N", "F"): -3, ("N", "P"): -2, ("N", "S"): 1, ("N", "T"): 0,
    ("N", "W"): -4, ("N", "Y"): -2, ("N", "V"): -3,
    ("D", "D"): 6, ("D", "C"): -3, ("D", "Q"): 0, ("D", "E"): 2, ("D", "G"): -1,
    ("D", "H"): -1, ("D", "I"): -3, ("D", "L"): -4, ("D", "K"): -1, ("D", "M"): -3,
    ("D", "F"): -3, ("D", "P"): -1, ("D", "S"): 0, ("D", "T"): -1, ("D", "W"): -4,
    ("D", "Y"): -3, ("D", "V"): -3,
    ("C", "C"): 9, ("C", "Q"): -3, ("C", "E"): -4, ("C", "G"): -3, ("C", "H"): -3,
    ("C", "I"): -1, ("C", "L"): -1, ("C", "K"): -3, ("C", "M"): -1, ("C", "F"): -2,
    ("C", "P"): -3, ("C", "S"): -1, ("C", "T"): -1, ("C", "W"): -2, ("C", "Y"): -2,
    ("C", "V"): -1,
    ("Q", "Q"): 5, ("Q", "E"): 2, ("Q", "G"): -2, ("Q", "H"): 0, ("Q", "I"): -3,
    ("Q", "L"): -2, ("Q", "K"): 1, ("Q", "M"): 0, ("Q", "F"): -3, ("Q", "P"): -1,
    ("Q", "S"): 0, ("Q", "T"): -1, ("Q", "W"): -2, ("Q", "Y"): -1, ("Q", "V"): -2,
    ("E", "E"): 5, ("E", "G"): -2, ("E", "H"): 0, ("E", "I"): -3, ("E", "L"): -3,
    ("E", "K"): 1, ("E", "M"): -2, ("E", "F"): -3, ("E", "P"): -1, ("E", "S"): 0,
    ("E", "T"): -1, ("E", "W"): -3, ("E", "Y"): -2, ("E", "V"): -2,
    ("G", "G"): 6, ("G", "H"): -2, ("G", "I"): -4, ("G", "L"): -4, ("G", "K"): -2,
    ("G", "M"): -3, ("G", "F"): -3, ("G", "P"): -2, ("G", "S"): 0, ("G", "T"): -2,
    ("G", "W"): -2, ("G", "Y"): -3, ("G", "V"): -3,
    ("H", "H"): 8, ("H", "I"): -3, ("H", "L"): -3, ("H", "K"): -1, ("H", "M"): -2,
    ("H", "F"): -1, ("H", "P"): -2, ("H", "S"): -1, ("H", "T"): -2, ("H", "W"): -2,
    ("H", "Y"): 2, ("H", "V"): -3,
    ("I", "I"): 4, ("I", "L"): 2, ("I", "K"): -3, ("I", "M"): 1, ("I", "F"): 0,
    ("I", "P"): -3, ("I", "S"): -2, ("I", "T"): -1, ("I", "W"): -3, ("I", "Y"): -1,
    ("I", "V"): 3,
    ("L", "L"): 4, ("L", "K"): -2, ("L", "M"): 2, ("L", "F"): 0, ("L", "P"): -3,
    ("L", "S"): -2, ("L", "T"): -1, ("L", "W"): -2, ("L", "Y"): -1, ("L", "V"): 1,
    ("K", "K"): 5, ("K", "M"): -1, ("K", "F"): -3, ("K", "P"): -1, ("K", "S"): 0,
    ("K", "T"): -1, ("K", "W"): -3, ("K", "Y"): -2, ("K", "V"): -2,
    ("M", "M"): 5, ("M", "F"): 0, ("M", "P"): -2, ("M", "S"): -1, ("M", "T"): -1,
    ("M", "W"): -1, ("M", "Y"): -1, ("M", "V"): 1,
    ("F", "F"): 6, ("F", "P"): -4, ("F", "S"): -2, ("F", "T"): -2, ("F", "W"): 1,
    ("F", "Y"): 3, ("F", "V"): -1,
    ("P", "P"): 7, ("P", "S"): -1, ("P", "T"): -1, ("P", "W"): -4, ("P", "Y"): -3,
    ("P", "V"): -2,
    ("S", "S"): 4, ("S", "T"): 1, ("S", "W"): -3, ("S", "Y"): -2, ("S", "V"): -2,
    ("T", "T"): 5, ("T", "W"): -2, ("T", "Y"): -2, ("T", "V"): 0,
    ("W", "W"): 11, ("W", "Y"): 2, ("W", "V"): -3,
    ("Y", "Y"): 7, ("Y", "V"): -1,
    ("V", "V"): 4,
}


def _blosum(a, b):
    a = a.upper()
    b = b.upper()
    if a == "X" or b == "X":
        return 0
    key = (a, b) if (a, b) in _BLOSUM62 else (b, a)
    return _BLOSUM62.get(key, -4)


def _parse_template(numbered_seq):
    positions = []
    for token in numbered_seq.split():
        match = re.match(r"^([A-Z])(\d+)([A-Z]?)$", token)
        if not match:
            raise ValueError("bad template token: %r" % token)
        aa, number, insertion = match.group(1), int(match.group(2)), match.group(3)
        positions.append((number, insertion, aa))
    return positions


# Kabat templates include explicit insertion slots (X + insertion code).
# Variable loops align into these slots instead of receiving ad-hoc numbers.
HEAVY_TEMPLATE = _parse_template(
    "E1 V2 Q3 L4 V5 E6 S7 G8 G9 G10 L11 V12 Q13 P14 G15 G16 S17 L18 R19 L20 "
    "S21 C22 A23 A24 S25 G26 F27 T28 F29 S30 "
    "S31 Y32 A33 M34 S35 X35A X35B W36 V37 R38 Q39 A40 P41 G42 K43 G44 L45 E46 W47 V48 A49 "
    "S50 I51 S52 X52A X52B X52C G53 G54 S55 T56 Y57 Y58 A59 D60 S61 V62 K63 G64 R65 "
    "F66 T67 I68 S69 R70 D71 N72 S73 K74 N75 T76 L77 Y78 L79 Q80 M81 N82 "
    "X82A X82B X82C S83 R84 A85 E86 D87 T88 A89 V90 Y91 "
    "C92 A93 K94 D95 Y96 G97 N98 Y99 Y100 "
    "X100A X100B X100C X100D X100E X100F X100G X100H X100I X100J X100K "
    "M101 D102 W103 G104 Q105 G106 T107 L108 V109 T110 V111 S112 S113"
)

KAPPA_TEMPLATE = _parse_template(
    "D1 I2 Q3 M4 T5 Q6 S7 P8 S9 S10 L11 S12 A13 S14 V15 G16 D17 R18 V19 T20 "
    "I21 T22 C23 R24 A25 S26 Q27 X27A X27B X27C X27D X27E X27F S28 V29 S30 S31 Y32 L33 A34 "
    "W35 Y36 Q37 Q38 K39 P40 G41 K42 A43 P44 K45 L46 L47 I48 Y49 "
    "A50 A51 S52 S53 L54 E55 S56 "
    "G57 V58 P59 S60 R61 F62 S63 G64 S65 G66 S67 G68 T69 D70 F71 T72 L73 "
    "T74 I75 S76 S77 L78 Q79 P80 E81 D82 F83 A84 T85 Y86 Y87 C88 "
    "Q89 Q90 Y91 N92 S93 L94 P95 X95A X95B X95C X95D X95E X95F Y96 T97 "
    "F98 G99 Q100 G101 T102 K103 V104 E105 I106 K107"
)

LAMBDA_TEMPLATE = _parse_template(
    "Q1 S2 V3 L4 T5 Q6 P7 P8 S9 A10 S11 G12 T13 P14 G15 Q16 R17 V18 T19 "
    "I20 T21 T22 C23 T24 G25 T26 S27 S28 X27A X27B X27C X27D X27E X27F D29 V30 "
    "G31 G32 Y33 N34 W35 Y36 Q37 Q38 H39 P40 G41 K42 A43 P44 K45 L46 M47 I48 Y49 "
    "D50 V51 S52 N53 R54 P55 S56 "
    "G57 V58 S59 N60 R61 F62 S63 G64 S65 K66 S67 G68 N69 T70 A71 S72 L73 "
    "T74 I75 S76 G77 L78 Q79 A80 E81 D82 E83 A84 D85 Y86 Y87 C88 "
    "S89 S90 Y91 T92 S93 S94 S95 X95A X95B X95C X95D X95E X95F T96 L97 "
    "F98 G99 G100 G101 T102 K103 L104 T105 V106 L107"
)

TEMPLATES = {
    "H": ("H", HEAVY_TEMPLATE),
    "K": ("L", KAPPA_TEMPLATE),
    "L": ("L", LAMBDA_TEMPLATE),
}


def _safe_name(value):
    value = re.sub(r"[^A-Za-z0-9_]", "_", str(value))
    return value.strip("_") or "cdr"


def _truthy(value):
    return str(value).strip().lower() not in ("0", "false", "no", "off", "")


def _split_chains(chains):
    text = str(chains).replace(",", "+").replace(" ", "")
    parts = [p for p in text.split("+") if p]
    if not parts:
        raise ValueError("chains is empty")
    return parts


def _looks_like_chains(value):
    """True if text looks like chain IDs (H+L, A, H/L) rather than an object name."""
    text = str(value).strip()
    if not text:
        return False
    text = text.replace(",", "+").replace(" ", "").replace("/", "+")
    parts = [p for p in text.split("+") if p]
    if not parts:
        return False
    return all(re.match(r"^[A-Za-z0-9]$", p) for p in parts)


def _enabled_objects():
    return list(cmd.get_object_list("enabled") or [])


def _resolve_object_and_chains(obj, chains):
    """Allow short forms: ab_cdr / ab_cdr antibody / ab_cdr H+L."""
    obj = str(obj or "").strip()
    chains = str(chains or "").strip()

    if obj and not chains and _looks_like_chains(obj):
        chains = obj.replace("/", "+")
        obj = ""

    if not obj:
        objects = _enabled_objects()
        if len(objects) == 0:
            raise ValueError("no enabled molecular object; load a structure first")
        if len(objects) > 1:
            raise ValueError(
                "multiple objects enabled; specify one, e.g. ab_cdr antibody"
            )
        obj = objects[0]
    return obj, chains


def _resi_token(resi, icode=""):
    icode = (icode or "").strip()
    if icode:
        return "%s%s" % (resi, icode)
    return str(resi)


def _nw_align(query, template, gap_open=-10, gap_extend=-1):
    """Semi-global NW: free end gaps on both sequences."""
    n = len(query)
    m = len(template)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]  # 0=diag, 1=up, 2=left

    for i in range(1, n + 1):
        score[i][0] = 0
        trace[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = 0
        trace[0][j] = 2

    for i in range(1, n + 1):
        qa = query[i - 1]
        for j in range(1, m + 1):
            ta = template[j - 1][2]
            diag = score[i - 1][j - 1] + _blosum(qa, ta)
            # affine-ish linear gap is enough here
            up = score[i - 1][j] + (gap_extend if trace[i - 1][j] == 1 else gap_open)
            left = score[i][j - 1] + (gap_extend if trace[i][j - 1] == 2 else gap_open)
            # Prefer match on ties to keep conserved anchors together.
            best = diag
            move = 0
            if up > best:
                best = up
                move = 1
            if left > best:
                best = left
                move = 2
            score[i][j] = best
            trace[i][j] = move

    # Best end: allow free terminal gaps.
    best_i, best_j, best_s = n, m, score[n][m]
    for i in range(n + 1):
        if score[i][m] > best_s:
            best_s = score[i][m]
            best_i, best_j = i, m
    for j in range(m + 1):
        if score[n][j] > best_s:
            best_s = score[n][j]
            best_i, best_j = n, j

    i, j = best_i, best_j
    pairs = []  # (query_index or None, template_index or None)
    while i > 0 and j > 0:
        move = trace[i][j]
        if move == 0:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif move == 1:
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    while i > 0:
        pairs.append((i - 1, None))
        i -= 1
    while j > 0:
        pairs.append((None, j - 1))
        j -= 1
    pairs.reverse()

    # Extend unmatched query tails beyond the best end as insertions.
    for qi in range(best_i, n):
        pairs.append((qi, None))

    return best_s, pairs


def _kabat_hcdr1_positions(length):
    if length <= 0:
        return []
    if length <= 5:
        return [(num, "") for num in range(31, 31 + length)]
    extra = length - 5
    if extra > 2:
        raise ValueError("HCDR1 too long for Kabat numbering (%d aa)" % length)
    return [(num, "") for num in range(31, 36)] + [
        (35, _KABAT_ALPHABET[i]) for i in range(extra)
    ]


def _kabat_hcdr2_positions(length):
    if length <= 0:
        return []
    if length <= 16:
        return [(num, "") for num in range(50, 50 + length)]
    extra = length - 16
    if extra > 3:
        raise ValueError("HCDR2 too long for Kabat numbering (%d aa)" % length)
    head = [(50, ""), (51, ""), (52, "")]
    inserts = [(52, _KABAT_ALPHABET[i]) for i in range(extra)]
    tail = [(num, "") for num in range(53, 66)]
    return (head + inserts + tail)[:length]


def _kabat_hcdr3_positions(length):
    if length <= 0:
        return []
    if length > 36:
        raise ValueError("HCDR3 too long for Kabat numbering (%d aa)" % length)
    insertions = max(length - 8, 0)
    ordered = [
        (100, ""), (99, ""), (98, ""), (97, ""), (96, ""), (95, ""),
        (101, ""), (102, ""),
    ]
    base = sorted(ordered[max(0, 8 - length):])
    return sorted(base + [(100, _KABAT_ALPHABET[i]) for i in range(insertions)])


def _kabat_lcdr1_positions(length):
    if length <= 0:
        return []
    insertions = max(length - 11, 0)
    if insertions > 6:
        raise ValueError("LCDR1 too long for Kabat numbering (%d aa)" % length)
    head = [(num, "") for num in range(24, 28)]
    mid = [(27, _KABAT_ALPHABET[i]) for i in range(insertions)]
    tail = [(num, "") for num in range(28, 35)]
    return (head + mid + tail)[:length]


def _kabat_lcdr2_positions(length):
    if length <= 0:
        return []
    if length > 7:
        raise ValueError("LCDR2 longer than Kabat definition (%d aa)" % length)
    return [(num, "") for num in range(50, 50 + length)]


def _kabat_lcdr3_positions(length):
    if length <= 0:
        return []
    if length > 35:
        raise ValueError("LCDR3 too long for Kabat numbering (%d aa)" % length)
    insertions = max(length - 9, 0)
    ordered = [
        (95, ""), (94, ""), (93, ""), (92, ""), (91, ""),
        (96, ""), (97, ""), (90, ""), (89, ""),
    ]
    base = sorted(ordered[max(0, 9 - length):])
    return sorted(base + [(95, _KABAT_ALPHABET[i]) for i in range(insertions)])


_KABAT_CDR_POSITION_FN = {
    "HCDR1": _kabat_hcdr1_positions,
    "HCDR2": _kabat_hcdr2_positions,
    "HCDR3": _kabat_hcdr3_positions,
    "LCDR1": _kabat_lcdr1_positions,
    "LCDR2": _kabat_lcdr2_positions,
    "LCDR3": _kabat_lcdr3_positions,
}

_KABAT_TEMPLATE_CDR_RANGES = {
    "H": {
        "HCDR1": (31, 35),
        "HCDR2": (50, 65),
        "HCDR3": (95, 102),
    },
    "L": {
        "LCDR1": (24, 34),
        "LCDR2": (50, 56),
        "LCDR3": (89, 97),
    },
}


def _cdr_name_for_template_num(chain_class, num, ins):
    ranges = _KABAT_TEMPLATE_CDR_RANGES[chain_class]
    for cdr_name, (lo, hi) in ranges.items():
        if lo <= num <= hi:
            return cdr_name
    return None


def _sequence_hcdr3_span(sequence):
    cys = [i for i, aa in enumerate(sequence) if aa == "C"]
    if len(cys) < 2:
        return None
    c2 = cys[1]
    for i in range(c2 + 3, len(sequence) - 1):
        if sequence[i] == "W" and sequence[i + 1] == "G":
            start = c2 + 3
            if start < i:
                return list(range(start, i))
    return None


def _sequence_lcdr3_span(sequence):
    cys = [i for i, aa in enumerate(sequence) if aa == "C"]
    if len(cys) < 2:
        return None
    c2 = cys[1]
    for i in range(c2 + 2, len(sequence) - 1):
        if sequence[i] == "F" and sequence[i + 1] == "G":
            start = c2 + 1
            if start < i:
                return list(range(start, i))
    return None


def _sequence_hcdr1_span(sequence):
    cys = [i for i, aa in enumerate(sequence) if aa == "C"]
    if not cys:
        return None
    c1 = cys[0]
    for i in range(c1 + 8, len(sequence) - 1):
        if sequence[i] == "W" and i + 1 < len(sequence) and sequence[i + 1] in "VI":
            start = max(c1 + 4, i - 5)
            if start < i:
                return list(range(start, i))
    return None


def _sequence_lcdr1_span(sequence):
    cys = [i for i, aa in enumerate(sequence) if aa == "C"]
    if not cys:
        return None
    c1 = cys[0]
    for i in range(c1 + 8, len(sequence) - 1):
        if sequence[i] == "W":
            end = i - 1
            while end > c1 + 1 and sequence[end] not in "VILMFYWA":
                end -= 1
            start = c1 + 1
            if start <= end:
                return list(range(start, end + 1))
    return None


def _collect_cdr_indices_from_pairs(pairs, template, chain_class, sequence=""):
    """Map CDR name -> ordered unique query indices using template alignment."""
    tagged = {}
    pair_by_qi = {}
    framework = {}
    for qi, tj in pairs:
        if qi is None:
            continue
        pair_by_qi[qi] = tj
        if tj is None:
            continue
        num, ins, _aa = template[tj]
        cdr_name = _cdr_name_for_template_num(chain_class, num, ins)
        if cdr_name is not None:
            tagged[qi] = cdr_name
        else:
            framework[qi] = num

    buckets = {name: [] for name in KABAT_CDR_RANGES[chain_class]}
    by_name = defaultdict(list)
    for qi, name in tagged.items():
        by_name[name].append(qi)

    cdr_limits = {
        "HCDR3": (95, 102, 103),
        "LCDR3": (89, 97, 98),
    }

    for cdr_name, hits in by_name.items():
        if not hits:
            continue
        lo, hi = min(hits), max(hits)
        limit = cdr_limits.get(cdr_name)
        if limit:
            start_num, end_num, fw_num = limit
            for qi, num in framework.items():
                if num >= fw_num:
                    hi = min(hi, qi - 1)
                if num < start_num and cdr_name.endswith("CDR3"):
                    lo = max(lo, qi + 1)
        selected = []
        for qi in range(lo, hi + 1):
            if qi in tagged and tagged[qi] == cdr_name:
                selected.append(qi)
            elif pair_by_qi.get(qi) is None:
                selected.append(qi)
        buckets[cdr_name] = selected

    if sequence:
        seq_spans = {
            "H": {
                "HCDR1": _sequence_hcdr1_span,
                "HCDR3": _sequence_hcdr3_span,
            },
            "L": {
                "LCDR1": _sequence_lcdr1_span,
                "LCDR3": _sequence_lcdr3_span,
            },
        }.get(chain_class, {})
        for cdr_name, fn in seq_spans.items():
            span = fn(sequence)
            if span:
                buckets[cdr_name] = span
    return buckets


def _apply_positions_to_indices(numbered, index_positions):
    pos_by_idx = dict(index_positions)
    out = []
    for idx, aa, num, ins in numbered:
        if idx in pos_by_idx:
            num, ins = pos_by_idx[idx]
            out.append((idx, aa, num, ins or ""))
        else:
            out.append((idx, aa, num, ins or ""))
    return out


def _fix_framework_after_cdr_renumber(numbered, cdr_indices, pairs, template, chain_class):
    """Drop stale CDR-range numbers from framework residues."""
    cdr_set = set()
    for indices in cdr_indices.values():
        cdr_set.update(indices)
    cdr_nums = set()
    for _name, (lo, hi) in KABAT_CDR_RANGES[chain_class].items():
        cdr_nums.update(range(lo, hi + 1))

    pair_num = {}
    for qi, tj in pairs:
        if qi is not None and tj is not None:
            num, ins, _aa = template[tj]
            pair_num[qi] = (num, ins or "")

    out = []
    for idx, aa, num, ins in numbered:
        ins = ins or ""
        if idx not in cdr_set and num in cdr_nums:
            if idx in pair_num:
                num, ins = pair_num[idx]
            else:
                num, ins = None, ""
        out.append((idx, aa, num, ins))
    return out
    pos_by_idx = dict(index_positions)
    out = []
    for idx, aa, num, ins in numbered:
        if idx in pos_by_idx:
            num, ins = pos_by_idx[idx]
            out.append((idx, aa, num, ins or ""))
        else:
            out.append((idx, aa, num, ins or ""))
    return out


def _renumber_kabat_cdrs(numbered, chain_class, cdr_indices):
    index_positions = {}
    for cdr_name, indices in cdr_indices.items():
        if not indices:
            continue
        fn = _KABAT_CDR_POSITION_FN.get(cdr_name)
        if fn is None:
            continue
        try:
            positions = fn(len(indices))
        except ValueError:
            continue
        for qi, pos in zip(indices, positions):
            index_positions[qi] = pos
    return _apply_positions_to_indices(numbered, index_positions)


def _kabat_insertion_site(chain_class, prev_num, next_num):
    if chain_class == "H":
        rules = [(31, 36, 35), (50, 66, 52), (66, 93, 82), (95, 103, 100)]
    else:
        rules = [(24, 35, 27), (50, 57, 52), (89, 98, 95)]
    for lo, hi, anchor in rules:
        if prev_num is not None and lo <= prev_num < hi:
            return anchor
        if next_num is not None and lo < next_num <= hi:
            return anchor
    return prev_num if prev_num is not None else next_num


def _assign_insertions(numbered, chain_class):
    """Assign Kabat insertion codes to unnumbered residues."""
    result = []
    pending = []
    last_num = None
    anchor_counts = defaultdict(int)

    def _flush_pending(next_num):
        nonlocal pending, last_num, anchor_counts
        if not pending:
            return
        anchor = _kabat_insertion_site(chain_class, last_num, next_num)
        if anchor is None:
            for pidx, paa in pending:
                result.append((pidx, paa, None, ""))
        else:
            for pidx, paa in pending:
                count = anchor_counts[anchor]
                if count >= len(_KABAT_ALPHABET):
                    result.append((pidx, paa, None, ""))
                    continue
                letter = _KABAT_ALPHABET[count]
                anchor_counts[anchor] += 1
                result.append((pidx, paa, anchor, letter))
        pending = []

    for idx, aa, num, ins in numbered:
        if num is None:
            pending.append((idx, aa))
            continue
        _flush_pending(num)
        result.append((idx, aa, num, ins or ""))
        last_num = num

    _flush_pending(None)
    return result


def _validate_ig_fv(sequence, chain_class, score, mapped):
    """Reject antigen/non-Ig chains using Fv sequence motifs, not exact Kabat slots."""
    seq = re.sub(r"[^A-Za-z]", "", sequence.upper())
    if mapped < 55:
        raise ValueError("alignment too weak (%d mapped positions)" % mapped)

    cys = [i for i, aa in enumerate(seq) if aa == "C"]
    if chain_class == "H":
        if len(cys) < 2:
            raise ValueError("not an Ig heavy chain (need 2 Cys)")
        if score < 80:
            raise ValueError("heavy-chain alignment score too low (%d)" % score)
        if not seq.startswith(("QVQL", "EVQL", "QMQL", "QLQL", "QVQ", "DVQL", "AVQL")):
            if score < 150:
                raise ValueError("sequence does not look like a heavy variable domain")
    else:
        if len(cys) < 2:
            raise ValueError("not an Ig light chain (need 2 Cys)")
        if score < 80:
            raise ValueError("light-chain alignment score too low (%d)" % score)
        if not seq.startswith(
            ("DIQMT", "EIVLT", "QSVLT", "SYELT", "DIVMT", "AIRMT", "QSALT")
        ):
            if score < 150:
                raise ValueError("sequence does not look like a light variable domain")


def _find_fv_slice(sequence, chain_class):
    """Return (start, end) slice covering the variable domain only."""
    seq = sequence.upper()
    n = len(seq)
    start = 0
    cys_positions = [i for i, aa in enumerate(seq) if aa == "C"]
    if not cys_positions:
        return 0, min(n, 120 if chain_class == "H" else 115)

    if chain_class == "H" and len(cys_positions) >= 2:
        c2 = cys_positions[1]
        fv_len = 120
    else:
        c2 = cys_positions[min(1, len(cys_positions) - 1)]
        fv_len = 115

    end = n
    for w in [i for i, aa in enumerate(seq) if aa == "W" and i > c2]:
        end = min(n, w + 11)
        break
    else:
        end = min(n, cys_positions[0] + fv_len)
    return start, end


def _trim_rows_to_fv(rows, chain_class):
    sequence = "".join(row["aa"] for row in rows)
    start, end = _find_fv_slice(sequence, chain_class)
    if start == 0 and end == len(rows):
        return rows, 0
    return rows[start:end], start


def _prefer_antibody_chains(annotated):
    """Keep antibody chains only when a complex also has antigen-like chains."""
    if len(annotated) <= 1:
        return annotated

    by_id = dict(annotated)
    preferred_ids = []
    if "H" in by_id:
        preferred_ids.append("H")
    for light_id in ("L", "K"):
        if light_id in by_id:
            preferred_ids.append(light_id)
            break
    if preferred_ids:
        return [(cid, by_id[cid]) for cid in preferred_ids]

    heavy = [(c, r) for c, r in annotated if r["chain_class"] == "H"]
    light = [(c, r) for c, r in annotated if r["chain_class"] == "L"]
    if heavy and light:
        heavy.sort(key=lambda x: x[1].get("score", 0), reverse=True)
        light.sort(key=lambda x: x[1].get("score", 0), reverse=True)
        return [heavy[0], light[0]]
    return annotated


def number_by_consensus(sequence):
    """Return (chain_class, scheme_label, score, numbered_residues, cdr_indices)."""
    sequence = "".join(AA1.get(x, x) if len(x) > 1 else x for x in sequence)
    sequence = re.sub(r"[^A-Za-z]", "", sequence.upper())
    if len(sequence) < 70:
        raise ValueError("sequence too short for a variable domain (%d aa)" % len(sequence))

    chain_class_hint = _guess_chain_class(sequence)
    start, end = _find_fv_slice(sequence, chain_class_hint)
    offset = start
    sequence = sequence[start:end]

    best = None
    for key, (chain_class, template) in TEMPLATES.items():
        score, pairs = _nw_align(sequence, template)
        numbered = []
        for qi, tj in pairs:
            if qi is None:
                continue
            if tj is None:
                numbered.append((qi, sequence[qi], None, ""))
            else:
                num, ins, _aa = template[tj]
                numbered.append((qi, sequence[qi], num, ins))
        numbered = _assign_insertions(numbered, chain_class)
        # Require conserved cysteines when possible.
        nums = {num for _i, _a, num, _ins in numbered if num is not None}
        bonus = 0
        if chain_class == "H":
            if 22 in nums:
                bonus += 50
            if 92 in nums:
                bonus += 50
            if 103 in nums:
                bonus += 30
        else:
            if 23 in nums:
                bonus += 50
            if 88 in nums:
                bonus += 50
        total = score + bonus
        item = (total, score, chain_class, key, numbered, pairs, template)
        if best is None or total > best[0]:
            best = item

    _total, score, chain_class, key, numbered, pairs, template = best
    cdr_indices = _collect_cdr_indices_from_pairs(pairs, template, chain_class, sequence)
    numbered = _renumber_kabat_cdrs(numbered, chain_class, cdr_indices)
    numbered = _fix_framework_after_cdr_renumber(
        numbered, cdr_indices, pairs, template, chain_class
    )
    if offset:
        numbered = [(qi + offset, aa, num, ins) for qi, aa, num, ins in numbered]
        cdr_indices = {
            name: [i + offset for i in idxs] for name, idxs in cdr_indices.items()
        }
    mapped = sum(1 for _i, _a, num, _ins in numbered if num is not None)
    if mapped < 60:
        raise ValueError(
            "consensus numbering failed (only %d positions mapped; score=%d)"
            % (mapped, score)
        )
    _validate_ig_fv(sequence, chain_class, score, mapped)
    return chain_class, "consensus:%s" % key, score, numbered, cdr_indices


def _guess_chain_class(sequence):
    sequence = sequence.upper()
    if sequence.startswith(("QVQL", "EVQL", "QMQL", "QLQL", "QVQ")):
        return "H"
    if sequence.startswith(("DIQMT", "EIVLT", "QSVLT", "SYELT", "QSVLT")):
        return "L"
    return "H"


def number_by_anarci(sequence, scheme="kabat"):
    """Return (chain_class, scheme_label, score, numbered_residues, cdr_indices)."""
    from anarci import anarci as run_anarci

    sequence = re.sub(r"[^A-Za-z]", "", sequence.upper())
    results = run_anarci(
        [("query", sequence)],
        scheme=scheme,
        output=False,
        allow_gap=True,
    )
    hit_list = results[0][0]
    if not hit_list:
        raise ValueError("ANARCI found no Ig domain")

    numbering, details = hit_list[0]
    start = int(details.get("query_start", 0))
    end = int(details.get("query_end", len(sequence) - 1))

    numbered = []
    for i in range(start):
        numbered.append((i, sequence[i], None, ""))

    seq_index = start
    for (num, ins), aa in numbering:
        if aa == "-":
            continue
        if seq_index > end:
            break
        numbered.append((seq_index, sequence[seq_index], int(num), (ins or "").strip()))
        seq_index += 1

    while seq_index < len(sequence):
        numbered.append((seq_index, sequence[seq_index], None, ""))
        seq_index += 1

    chain_type = details.get("chain_type", "H")
    chain_class = "H" if chain_type == "H" else "L"
    cdr_indices = {}
    for cdr_name in cdr_names_for_class(chain_class):
        lo, hi = KABAT_CDR_RANGES[chain_class][cdr_name]
        cdr_indices[cdr_name] = [
            idx for idx, _aa, num, _ins in numbered if num is not None and lo <= num <= hi
        ]
    score = float(details.get("bit_score") or 0)
    mapped = sum(1 for _i, _a, num, _ins in numbered if num is not None)
    fv_seq = "".join(aa for _i, aa, num, _ins in numbered if num is not None)
    if fv_seq:
        _validate_ig_fv(fv_seq, chain_class, max(score, 200), mapped)
    return chain_class, "anarci:%s" % scheme, score, numbered, cdr_indices


def number_sequence(sequence, scheme="kabat", engine="auto"):
    engine = str(engine).strip().lower()
    errors = []
    if engine in ("auto", "anarci"):
        try:
            return number_by_anarci(sequence, scheme=scheme)
        except Exception as exc:
            errors.append("anarci: %s" % exc)
            if engine == "anarci":
                raise ValueError(
                    "ANARCI numbering failed (%s). Install HMMER and ensure "
                    "`hmmscan` is on PATH, or use engine=consensus."
                    % exc
                )
    if engine in ("auto", "consensus"):
        if scheme.lower() != "kabat":
            raise ValueError("built-in consensus engine currently supports only Kabat")
        return number_by_consensus(sequence)
    raise ValueError("unknown engine %r (use auto, anarci, or consensus)" % engine)


def cdr_names_for_class(chain_class):
    return list(KABAT_CDR_RANGES[chain_class].keys())


def residues_in_cdr(numbered, chain_class):
    """Map CDR name -> list of seq indices."""
    buckets = {name: [] for name in cdr_names_for_class(chain_class)}
    ranges = KABAT_CDR_RANGES[chain_class]
    for idx, _aa, num, _ins in numbered:
        if num is None:
            continue
        for name, (lo, hi) in ranges.items():
            if lo <= num <= hi:
                buckets[name].append(idx)
                break
    return buckets


def extract_chain_residues(selection, state=1):
    """Return ordered list of dicts for polymer residues in selection."""
    if cmd is None:
        raise RuntimeError("PyMOL cmd is not available")

    raw = {"stored_cdr_atoms": []}
    namespace = {"stored_cdr_atoms": raw["stored_cdr_atoms"]}
    cmd.iterate_state(
        state,
        "(%s) and polymer.protein and name CA" % selection,
        "stored_cdr_atoms.append((segi, chain, resi, resn, resv))",
        space=namespace,
    )
    rows = []
    seen = set()
    for segi, chain, resi, resn, resv in raw["stored_cdr_atoms"]:
        key = (segi, chain, resi)
        if key in seen:
            continue
        seen.add(key)
        aa = AA1.get(resn, "X")
        # resi may include insertion code as trailing letter
        match = re.match(r"^(-?\d+)([A-Za-z]?)$", str(resi))
        if match:
            resi_num, icode = match.group(1), match.group(2)
        else:
            resi_num, icode = str(resi), ""
        rows.append(
            {
                "segi": segi or "",
                "chain": chain or "",
                "resi": str(resi),
                "resi_num": resi_num,
                "icode": icode,
                "resn": resn,
                "aa": aa,
                "resv": resv,
            }
        )
    if not rows:
        raise ValueError("no CA atoms found in selection: %s" % selection)
    return rows


def _residue_expr(obj_sel, residue):
    parts = ["(%s)" % obj_sel, "resi %s" % residue["resi"]]
    if residue["segi"]:
        parts.append("segi %s" % residue["segi"])
    if residue["chain"]:
        parts.append("chain %s" % residue["chain"])
    else:
        parts.append('chain ""')
    return "(%s)" % " and ".join(parts)


def _build_selection_expression(obj_sel, residues):
    if not residues:
        return "none"
    # Group by chain/segi for compact expressions.
    groups = defaultdict(list)
    for res in residues:
        groups[(res["segi"], res["chain"])].append(res["resi"])
    chunks = []
    for (segi, chain), resis in groups.items():
        parts = ["(%s)" % obj_sel]
        if segi:
            parts.append("segi %s" % segi)
        if chain:
            parts.append("chain %s" % chain)
        else:
            parts.append('chain ""')
        parts.append("resi " + "+".join(resis))
        chunks.append("(%s)" % " and ".join(parts))
    return " or ".join(chunks)


def annotate_chain_rows(rows, scheme="kabat", engine="auto"):
    chain_class_hint = _guess_chain_class("".join(row["aa"] for row in rows))
    trimmed_rows, _offset = _trim_rows_to_fv(rows, chain_class_hint)
    sequence = "".join(row["aa"] for row in trimmed_rows)
    chain_class, engine_label, score, numbered, cdr_indices = number_sequence(
        sequence, scheme=scheme, engine=engine
    )
    index_to_row = {i: trimmed_rows[i] for i in range(len(trimmed_rows))}
    cdr_residues = {}
    for name, indices in cdr_indices.items():
        cdr_residues[name] = [index_to_row[i] for i in indices if i in index_to_row]
    return {
        "chain_class": chain_class,
        "engine": engine_label,
        "score": score,
        "sequence": sequence,
        "numbered": numbered,
        "cdr_indices": cdr_indices,
        "cdr_residues": cdr_residues,
        "fv_rows": trimmed_rows,
    }


def ab_cdr(
    obj="",
    chains="",
    scheme="kabat",
    engine="consensus",
    prefix="cdr",
    state=1,
    sticks=0,
    quiet=0,
):
    """
DESCRIPTION

    Color Kabat CDR1/2/3 only. Other residues keep their current display.

USAGE

    ab_cdr
    ab_cdr antibody
    ab_cdr H+L
    ab_cdr antibody, H+L

ARGUMENTS

    object
        Molecular object. Empty = the single enabled object.
        If this looks like chain IDs (H+L / A), it is treated as chains.

    chains
        Chain IDs, e.g. H+L or A. Empty = try all protein chains and
        skip non-Ig chains.

    sticks
        Also show CDR sticks (1/0). Default: 0 (color only).
    """
    if cmd is None:
        raise RuntimeError("PyMOL is required for ab_cdr")

    quiet = _truthy(quiet)
    sticks = _truthy(sticks)
    scheme = str(scheme).strip().lower() or "kabat"
    if scheme != "kabat":
        raise ValueError("this version currently supports scheme=kabat only")

    obj, chains = _resolve_object_and_chains(obj, chains)
    prefix = _safe_name(prefix)
    state = int(state)

    if chains:
        chain_ids = _split_chains(chains)
    else:
        chain_ids = []
        space = {"cset": set()}
        cmd.iterate(
            "(%s) and polymer.protein and name CA" % obj,
            "cset.add(chain)",
            space=space,
        )
        chain_ids = sorted([c for c in space["cset"] if c is not None and c != ""])
        if not chain_ids:
            chain_ids = [""]

    existing = cmd.get_names("all")
    for name in existing:
        if name == prefix or name.startswith(prefix + "_"):
            cmd.delete(name)

    base_sel = "(%s) and polymer.protein" % obj
    cmd.select(prefix + "_fv", "none")
    all_cdr_expr_parts = []
    summaries = []
    union_created = set()
    annotated = []
    skipped = []

    for chain_id in chain_ids:
        if chain_id:
            chain_sel = "(%s) and chain %s" % (base_sel, chain_id)
        else:
            chain_sel = '(%s) and chain ""' % base_sel
        if cmd.count_atoms(chain_sel) == 0:
            if not quiet:
                print("skip empty chain %r" % chain_id)
            continue
        rows = extract_chain_residues(chain_sel, state=state)
        try:
            result = annotate_chain_rows(rows, scheme=scheme, engine=engine)
        except Exception as exc:
            skipped.append((chain_id or "?", str(exc)))
            continue
        annotated.append((chain_id or "X", result))

    # Complexes: keep antibody/nanobody chains only; drop antigen-like hits.
    if not chains:
        annotated = _prefer_antibody_chains(annotated)

    if not annotated:
        if skipped:
            detail = "; ".join("%s: %s" % (c, e) for c, e in skipped)
            raise ValueError(
                "no antibody/nanobody variable domains were annotated; %s"
                % detail
            )
        raise ValueError(
            "no antibody/nanobody variable domains were annotated; "
            "check chain IDs or try engine=consensus"
        )

    for chain_tag, result in annotated:
        fv_sel = _build_selection_expression(obj, result.get("fv_rows") or [])
        if fv_sel != "none":
            cmd.select(
                prefix + "_fv",
                "(%s) or (%s)" % (prefix + "_fv", fv_sel),
            )

        for cdr_name, residues in result["cdr_residues"].items():
            sel_name = "%s_%s_%s" % (prefix, chain_tag, cdr_name)
            expr = _build_selection_expression(obj, residues)
            cmd.select(sel_name, expr if residues else "none")
            if residues:
                all_cdr_expr_parts.append(sel_name)
            summaries.append(
                (chain_tag, result["chain_class"], cdr_name, len(residues), result["engine"])
            )

        for cdr_name in result["cdr_residues"]:
            union_name = "%s_%s" % (prefix, cdr_name)
            chain_sel_name = "%s_%s_%s" % (prefix, chain_tag, cdr_name)
            if union_name in union_created:
                cmd.select(union_name, "(%s) or (%s)" % (union_name, chain_sel_name))
            else:
                cmd.select(union_name, chain_sel_name)
                union_created.add(union_name)

    if all_cdr_expr_parts:
        cmd.select(prefix + "_CDRs", " or ".join("(%s)" % n for n in all_cdr_expr_parts))
    else:
        cmd.select(prefix + "_CDRs", "none")

    # Color CDRs only; leave the rest of the scene unchanged.
    for cdr_name, color in CDR_COLORS.items():
        sel = "%s_%s" % (prefix, cdr_name)
        if sel in cmd.get_names("all") and cmd.count_atoms(sel) > 0:
            cmd.color(color, sel)
            if sticks:
                cmd.show("sticks", sel)

    if not quiet:
        print("ab_cdr %s  scheme=%s  object=%s" % (__version__, scheme, obj))
        kept = "+".join(sorted({c for c, _r in annotated}))
        print("antibody chains: %s" % kept)
        for chain_tag, chain_class, cdr_name, count, eng in summaries:
            print(
                "  chain %s (%s) %-5s  %2d residues  [%s]"
                % (chain_tag, chain_class, cdr_name, count, eng)
            )
        if skipped and not chains:
            print(
                "skipped non-Ig/antigen chains: %s"
                % "; ".join("%s: %s" % (c, e) for c, e in skipped)
            )

    return summaries


def ab_cdr_help():
    """Print short usage help."""
    print(
        """
ab_cdr / cdr — color Kabat CDRs only (rest unchanged)

  ab_cdr              # one enabled object; auto antibody chains
  ab_cdr complex      # antigen–antibody complex: antibody only
  ab_cdr H+L          # force specific chains
  ab_cdr antibody, H+L
  cdr                 # same as ab_cdr

Notes:
  - Only CDR residues are recolored. Framework/display stay as-is.
  - Complexes: antigen chains are skipped automatically when chains are omitted.
  - Default sticks=0. Use sticks=1 if you also want CDR sticks.
  - Fab: H+L (or omit). Nanobody: single VHH chain (or omit).
""".strip()
    )


if cmd is not None:
    cmd.extend("ab_cdr", ab_cdr)
    cmd.extend("cdr", ab_cdr)
    cmd.extend("ab_cdr_help", ab_cdr_help)
    print("Loaded commands: ab_cdr, cdr, ab_cdr_help")
