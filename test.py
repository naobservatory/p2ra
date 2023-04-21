#!/usr/bin/env python3

import unittest

import pathogens
from pathogen_properties import *


class TestPathogens(unittest.TestCase):
    def test_hsv1_imported(self):
        self.assertIn("hsv_1", pathogens.pathogens)

    def test_properties_exist(self):
        for pathogen_name, pathogen in pathogens.pathogens.items():
            with self.subTest(pathogen=pathogen_name):
                self.assertIsInstance(pathogen.background, str)
                self.assertIsInstance(pathogen.pathogen_chars, PathogenChars)

                for estimate in pathogen.estimate_prevalences():
                    self.assertIsInstance(estimate, Prevalence)
>>>>>>> jefftk-normal-vars


if __name__ == "__main__":
    unittest.main()
