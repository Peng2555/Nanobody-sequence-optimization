# -*- coding: utf-8 -*-
"""Unit tests for Kabat CDR consensus numbering (no PyMOL required)."""

from __future__ import print_function

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import pymol_ab_cdr_annotator as cdr


VH = (
    "QVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGR"
    "FTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS"
)
VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQSVSSAVAWYQQKPGKAPKLLIYSASSLYSGVPSRFSGSRS"
    "GTDFTLTISSLQPEDFATYYCQQYNSLPYTFGQGTKVEIK"
)
VHH = (
    "QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGR"
    "FTISRDNAKNTVYLQMNSLKPEDTAVYYCAADSSGSYYYTGGSYYWYWGQGTQVTVSS"
)
VL_LAMBDA = (
    "QSVLTQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNSNRPSGVPDRFSGSK"
    "SGTSASLAITGLQAEDEADYYCQSYDSSLSGYVFGTGTKVTVL"
)


def _label(numbered, idx):
    for i, aa, num, ins in numbered:
        if i == idx:
            return "%s%s%s" % (aa, num, ins)
    raise KeyError(idx)


def _cdr_seq(numbered, cdr_indices, cdr_name):
    seq = dict((i, aa) for i, aa, num, ins in numbered)
    return "".join(seq[i] for i in cdr_indices.get(cdr_name, []))


class TestKabatCdrConsensus(unittest.TestCase):
    def test_looks_like_chains(self):
        self.assertTrue(cdr._looks_like_chains("H+L"))
        self.assertTrue(cdr._looks_like_chains("A"))
        self.assertTrue(cdr._looks_like_chains("H/L"))
        self.assertFalse(cdr._looks_like_chains("antibody"))
        self.assertFalse(cdr._looks_like_chains("1mlc"))

    def test_prefer_hl_over_antigen_chain_ids(self):
        annotated = [
            ("H", {"chain_class": "H", "score": 400}),
            ("L", {"chain_class": "L", "score": 400}),
            ("A", {"chain_class": "H", "score": 130}),
        ]
        kept = cdr._prefer_antibody_chains(annotated)
        self.assertEqual([c for c, _ in kept], ["H", "L"])

    def test_non_ig_sequence_rejected(self):
        # Lysozyme-like stretch: should not pass Ig Cys/Trp checks.
        decoy = (
            "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
        )
        with self.assertRaises(ValueError):
            cdr.number_sequence(decoy, engine="consensus")

    def test_templates_have_unique_positions(self):
        for name, (_cls, template) in cdr.TEMPLATES.items():
            keys = [(num, ins) for num, ins, _aa in template]
            self.assertEqual(len(keys), len(set(keys)), name)

    def test_heavy_cdrs(self):
        chain_class, _engine, _score, numbered, cdr_indices = cdr.number_sequence(
            VH, engine="consensus"
        )
        self.assertEqual(chain_class, "H")
        self.assertEqual(_cdr_seq(numbered, cdr_indices, "HCDR1"), "SYAMS")
        self.assertEqual(
            _cdr_seq(numbered, cdr_indices, "HCDR2"),
            "AISGSGGSTYYADSVKGR",
        )
        h3 = _cdr_seq(numbered, cdr_indices, "HCDR3")
        self.assertTrue(h3.startswith("DIQY"))
        self.assertIn("MDV", h3)
        self.assertEqual(len(h3), 13)
        labels = dict((i, (n, ins)) for i, aa, n, ins in numbered)
        self.assertNotIn(95, [n for n, _ins in labels.values() if _ins in ("A", "B", "C")])
        self.assertEqual(labels[110], (102, ""))
        self.assertTrue(any(n == 100 and ins for n, ins in labels.values()))

    def test_light_cdrs(self):
        chain_class, _engine, _score, numbered, cdr_indices = cdr.number_sequence(
            VL, engine="consensus"
        )
        self.assertEqual(chain_class, "L")
        self.assertEqual(_cdr_seq(numbered, cdr_indices, "LCDR1"), "RASQSVSSAVA")
        self.assertEqual(_cdr_seq(numbered, cdr_indices, "LCDR2"), "SASSLYS")
        self.assertEqual(_cdr_seq(numbered, cdr_indices, "LCDR3"), "QQYNSLPYT")

    def test_lambda_light_template_has_position_22(self):
        nums = [num for num, ins, aa in cdr.LAMBDA_TEMPLATE]
        self.assertIn(22, nums)

    def test_lambda_light_chain(self):
        chain_class, engine, _score, numbered, cdr_indices = cdr.number_sequence(
            VL_LAMBDA, engine="consensus"
        )
        self.assertEqual(chain_class, "L")
        self.assertTrue(engine.endswith(":L"))
        self.assertEqual(_cdr_seq(numbered, cdr_indices, "LCDR1"), "TGSSSNIGAGYDV")
        self.assertEqual(_cdr_seq(numbered, cdr_indices, "LCDR3"), "QSYDSSLSGYV")

    def test_nanobody_only_heavy_cdrs(self):
        chain_class, _engine, _score, numbered, cdr_indices = cdr.number_sequence(
            VHH, engine="consensus"
        )
        self.assertEqual(chain_class, "H")
        self.assertEqual(set(cdr_indices), {"HCDR1", "HCDR2", "HCDR3"})
        self.assertGreaterEqual(len(cdr_indices["HCDR3"]), 10)
        self.assertEqual(_cdr_seq(numbered, cdr_indices, "HCDR1"), "SYAMG")

    def test_constant_domain_tail_still_maps_anchors(self):
        long_vh = VH + "ASTKGPSVFPLAPSSKSTSGGTAALGCLVK"
        _cls, _engine, _score, numbered, _cdr_indices = cdr.number_sequence(
            long_vh, engine="consensus"
        )
        by_num = {(num, ins): aa for _i, aa, num, ins in numbered if num is not None}
        self.assertEqual(by_num.get((22, "")), "C")
        self.assertEqual(by_num.get((92, "")), "C")
        self.assertEqual(by_num.get((103, "")), "W")

    def test_fv_trim_excludes_constant_domain_from_numbering(self):
        long_vh = VH + "ASTKGPSVFPLAPSSKSTSGGTAALGCLVK"
        rows = [{"aa": aa, "resi": str(i + 1), "chain": "H", "segi": "", "icode": ""} for i, aa in enumerate(long_vh)]
        result = cdr.annotate_chain_rows(rows, engine="consensus")
        self.assertLess(len(result["fv_rows"]), len(rows))
        self.assertEqual(
            "".join(r["aa"] for r in result["fv_rows"]),
            VH,
        )


if __name__ == "__main__":
    unittest.main()
