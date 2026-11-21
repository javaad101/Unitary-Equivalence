import unittest
from bell_ricochet import (
    hilbert_schmidt_overlap,
    l2_norm_diff,
    postselected_bell_overlap,
    random_unitary,
    ricochet_lhs,
    ricochet_rhs,
    run_contextual_equivalence_demo,
)


class BellRicochetTests(unittest.TestCase):
    def test_ricochet_property_holds_to_machine_precision(self):
        for d in (2, 4, 8, 12):
            A = random_unitary(d, seed=100 + d)
            self.assertLess(l2_norm_diff(ricochet_lhs(A), ricochet_rhs(A)), 1e-13)

    def test_postselected_bell_equals_hilbert_schmidt(self):
        U1 = random_unitary(8, seed=1)
        U2 = random_unitary(8, seed=2)
        same = postselected_bell_overlap(U1, U1)
        self.assertAlmostEqual(same.real, 1.0, places=12)
        self.assertAlmostEqual(same.imag, 0.0, places=12)

        pair = postselected_bell_overlap(U1, U2)
        expected = hilbert_schmidt_overlap(U1, U2)
        self.assertAlmostEqual(pair.real, expected.real, places=12)
        self.assertAlmostEqual(pair.imag, expected.imag, places=12)

    def test_unitary_preserves_distinction_erased_by_scalar_projection(self):
        res = run_contextual_equivalence_demo()
        self.assertEqual(res["scalar_projection_a"], res["scalar_projection_b"])
        self.assertEqual(res["ctx1_yield_a"], res["ctx1_yield_b"])
        self.assertEqual(res["ctx2_yield_a"], 1.0)
        self.assertEqual(res["ctx2_yield_b"], 0.0)
        self.assertAlmostEqual(
            res["l2_distance_before_unitary"],
            res["l2_distance_after_unitary"],
            places=12,
        )


if __name__ == "__main__":
    unittest.main()
