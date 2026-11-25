"""Verify geometry matching and deletion guards with small numeric fixtures."""

import unittest

from delivery import select_occurrences, pair_inventory, pair_subset


def solid(x):
    return {"tight_box":[x,0,0,x+3,4,5],"center":[x+1.5,2,2.5],"volume":60,"faces":6}


class ExportBoundary(unittest.TestCase):
    def test_unique_current_target_uses_geometry(self):
        remove, review = select_occurrences([solid(17),solid(6)], {"I00001":[solid(6)]})
        self.assertEqual(remove,{"I00001":[2]})
        self.assertFalse(review)

    def test_absent_and_coincident_geometry_are_held(self):
        remove, review = select_occurrences([solid(6),solid(6)], {"I00001":[solid(6)],"I00002":[solid(8)]})
        self.assertFalse(remove)
        self.assertEqual(set(review),{"I00001","I00002"})

    def test_overlap_between_candidate_actions_is_held(self):
        remove, review = select_occurrences([solid(6)], {"I00001":[solid(6)],"I00002":[solid(6)]})
        self.assertFalse(remove)
        self.assertEqual(set(review),{"I00001","I00002"})

    def test_protected_body_holds_entire_multibody_candidate(self):
        remove, review = select_occurrences([solid(6),solid(8)], {"I00001":[solid(6),solid(8)]}, {2})
        self.assertFalse(remove)
        self.assertIn("I00001",review)

    def test_split_surface_needs_independent_area_agreement(self):
        target = dict(solid(6),faces=8,area=94.)
        same_surface = dict(solid(6),faces=7,area=94.)
        other_surface = dict(same_surface,area=95.)
        self.assertEqual(select_occurrences([same_surface],{"I00001":[target]})[0],{"I00001":[1]})
        self.assertFalse(select_occurrences([other_surface],{"I00001":[target]})[0])
        with self.assertRaises(ValueError):
            pair_inventory([same_surface],[target])

    def test_equal_count_can_still_be_wrong(self):
        with self.assertRaises(ValueError):
            pair_inventory([solid(1),solid(2)],[solid(2),solid(3)])

    def test_inventory_preserves_multiplicity_and_ignores_order(self):
        proof = pair_inventory([solid(1),solid(1),solid(8)],[solid(8),solid(1),solid(1)])
        self.assertEqual(set(proof.values()),{1,2,3})
        self.assertEqual(proof[3],1)

    def test_subset_can_verify_repeated_protected_body(self):
        proof = pair_subset([solid(1),solid(1)],[solid(8),solid(1),solid(9),solid(1)])
        self.assertEqual(set(proof.values()),{2,4})


if __name__ == "__main__":
    unittest.main()
