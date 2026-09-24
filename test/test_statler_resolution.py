import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"script"))
from analyze_statler_resolution import rebound_metrics, reference_force_bias


class ResolutionTest(unittest.TestCase):
    def test_reference_correction_refines(self):
        coarse = np.array(reference_force_bias(500)["correction_over_initial_gravity"])
        fine = np.array(reference_force_bias(1000)["correction_over_initial_gravity"])
        self.assertLess(coarse[1], 0)
        self.assertGreater(coarse[3], 0)
        self.assertTrue(np.all(np.abs(fine) < np.abs(coarse)))
        self.assertAlmostEqual(coarse[1], -.09441873, places=7)

    def test_monotone_and_rebound(self):
        t = np.geomspace(1, 10000, 10000)
        decreasing = t**(-.2)
        result = rebound_metrics(t, decreasing)
        self.assertEqual(result["positive_log_variation"], 0)
        self.assertEqual(result["rebound_times_trh"], [])
        bumped = decreasing*np.exp(.8*np.exp(-((np.log(t)-np.log(300))/.2)**2))
        result = rebound_metrics(t, bumped)
        self.assertEqual(len(result["rebound_times_trh"]), 1)
        self.assertGreater(result["rebound_prominence_fraction"][0], .5)

    def test_invalid_histories(self):
        with self.assertRaises(ValueError):
            rebound_metrics([0, 100, 99], [1, 2, 3])
        with self.assertRaises(ValueError):
            rebound_metrics([0, 100], [1, 0])
        with self.assertRaises(ValueError):
            rebound_metrics([0, 10], [1, 2])


if __name__ == "__main__":
    unittest.main()
