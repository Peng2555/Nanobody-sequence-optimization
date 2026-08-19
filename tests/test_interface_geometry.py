import unittest

import pymol_ab_interface_analyzer as analyzer


def atom(index, x, y, z, resn="ALA", resi="1", name="CB", chain="H", elem="C"):
    return {
        "model": "complex",
        "index": index,
        "segi": "",
        "chain": chain,
        "resi": resi,
        "resn": resn,
        "name": name,
        "alt": "",
        "elem": elem,
        "coord": (float(x), float(y), float(z)),
        "vdw": analyzer.ELEMENT_VDW.get(elem, 1.7),
        "b": 0.0,
    }


class GeometryTests(unittest.TestCase):
    def test_distance_and_centroid(self):
        self.assertAlmostEqual(analyzer._distance_xyz((0, 0, 0), (3, 4, 0)), 5.0)
        self.assertEqual(analyzer._centroid([(0, 0, 0), (2, 4, 6)]), (1.0, 2.0, 3.0))

    def test_acute_angle_is_direction_invariant(self):
        self.assertAlmostEqual(analyzer._acute_angle_degrees((0, 0, 1), (0, 0, -2)), 0.0)
        self.assertAlmostEqual(analyzer._acute_angle_degrees((1, 0, 0), (0, 1, 0)), 90.0)

    def test_grid_pairs_obeys_cutoff(self):
        left = [atom(1, 0, 0, 0), atom(2, 20, 20, 20, resi="2")]
        right = [
            atom(3, 3.99, 0, 0, chain="A"),
            atom(4, 4.01, 0, 0, resi="2", chain="A"),
        ]
        pairs = analyzer._grid_pairs(left, right, 4.0)
        self.assertEqual(len(pairs), 1)
        self.assertAlmostEqual(pairs[0][2], 3.99)

    def test_parallel_aromatic_geometry(self):
        store = {}
        ring1 = {
            "residue": ("", "H", "10", "TYR"),
            "centroid": (0.0, 0.0, 0.0),
            "normal": (0.0, 0.0, 1.0),
            "atom": atom(1, 0, 0, 0, "TYR", "10", "CG", "H"),
        }
        ring2 = {
            "residue": ("", "A", "20", "PHE"),
            "centroid": (0.0, 0.0, 4.5),
            "normal": (0.0, 0.0, -1.0),
            "atom": atom(2, 0, 0, 4.5, "PHE", "20", "CG", "A"),
        }
        analyzer._ring_interactions([ring1], [ring2], store, analyzer.DEFAULTS)
        self.assertIn(("pi_stacking", ring1["residue"], ring2["residue"]), store)

    def test_t_shaped_aromatic_geometry(self):
        store = {}
        ring1 = {
            "residue": ("", "H", "10", "TYR"),
            "centroid": (0.0, 0.0, 0.0),
            "normal": (0.0, 0.0, 1.0),
            "atom": atom(1, 0, 0, 0, "TYR", "10", "CG", "H"),
        }
        ring2 = {
            "residue": ("", "A", "20", "PHE"),
            "centroid": (0.0, 0.0, 4.0),
            "normal": (1.0, 0.0, 0.0),
            "atom": atom(2, 0, 0, 4.0, "PHE", "20", "CG", "A"),
        }
        analyzer._ring_interactions([ring1], [ring2], store, analyzer.DEFAULTS)
        self.assertIn(("t_stacking", ring1["residue"], ring2["residue"]), store)

    def test_cation_pi_geometry(self):
        store = {}
        cation_atom = atom(1, 0, 0, 4.0, "LYS", "10", "NZ", "H", "N")
        cation = {
            "residue": analyzer._residue_id(cation_atom),
            "center": cation_atom["coord"],
            "atom": cation_atom,
        }
        ring_atom = atom(2, 0, 0, 0, "TYR", "20", "CG", "A")
        ring = {
            "residue": analyzer._residue_id(ring_atom),
            "centroid": (0.0, 0.0, 0.0),
            "normal": (0.0, 0.0, 1.0),
            "atom": ring_atom,
        }
        analyzer._cation_pi_interactions([cation], [ring], True, store, analyzer.DEFAULTS)
        self.assertIn(("cation_pi", cation["residue"], ring["residue"]), store)


if __name__ == "__main__":
    unittest.main()
