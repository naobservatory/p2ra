#!/usr/bin/env python3

import unittest

import pathogens
from pathogen_properties import *

class TestPathogens(unittest.TestCase):
    def test_hsv1_imported(self):
        self.assertIn("hsv_1",  pathogens.pathogens)

    def test_properties_exist(self):
        for pathogen in pathogens.pathogens:
            with self.subTest(pathogen=pathogen):
                self.assertIsInstance(
                    pathogens.pathogens[pathogen].background, str)

                self.assertIsInstance(
                    pathogens.pathogens[pathogen].pathogen_chars,
                    PathogenChars)

                for prevalence_var_name, prevalence_var in pathogens.pathogens[
                        pathogen].prevalence_vars.items():
                    self.assertIsInstance(prevalence_var, PrevalenceVariable)


if __name__ == "__main__":
    unittest.main()

