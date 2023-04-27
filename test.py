#!/usr/bin/env python3

import unittest
from collections import Counter

import mgs
import pathogens
from mgs import MGSData
from pathogen_properties import *
from tree import Tree


class TestPathogens(unittest.TestCase):
    def test_hsv1_imported(self):
        self.assertIn("hsv_1", pathogens.pathogens)

    def test_summarize_location(self):
        us_2019, la_2020 = pathogens.pathogens["hiv"].estimate_prevalences()
        self.assertEqual(us_2019.summarize_location(), "United States")
        self.assertEqual(us_2019.summarize_date(), "2019")

        self.assertEqual(
            la_2020.summarize_location(),
            "Los Angeles, California, United States",
        )
        self.assertEqual(la_2020.summarize_date(), "2020")

    def test_properties_exist(self):
        for pathogen_name, pathogen in pathogens.pathogens.items():
            with self.subTest(pathogen=pathogen_name):
                self.assertIsInstance(pathogen.background, str)

                self.assertIsInstance(pathogen.pathogen_chars, PathogenChars)

                for estimate in pathogen.estimate_prevalences():
                    self.assertIsInstance(estimate, Prevalence)


class TestMGS(unittest.TestCase):
    repo = mgs.GitHubRepo(
        user="naobservatory",
        repo="mgs-pipeline",
        branch="47e2025f35168d3f414ae62928f6a14dd3f7c23d",
    )

    def test_load_bioprojects(self):
        bps = mgs.load_bioprojects(self.repo)
        # Rothman
        self.assertIn(mgs.BioProject("PRJNA729801"), bps)

    def test_load_sample_attributes(self):
        samples = mgs.load_sample_attributes(self.repo)
        # Randomly picked Rothman sample
        self.assertIn(mgs.Sample("SRR14530726"), samples)

    def test_load_sample_counts(self):
        sample_counts = mgs.load_sample_counts(self.repo)
        for p in ["sars_cov_2", "hiv"]:
            with self.subTest(pathogen=p):
                self.assertIn(
                    pathogens.pathogens[p].pathogen_chars.taxid,
                    sample_counts,
                )

    def test_load_tax_tree(self):
        tree = mgs.load_tax_tree(self.repo)
        for p in ["sars_cov_2", "hiv", "norovirus"]:
            with self.subTest(pathogen=p):
                self.assertIn(
                    pathogens.pathogens[p].pathogen_chars.taxid, tree
                )

    def test_count_reads(self):
        taxtree = Tree(mgs.TaxID(0), [Tree(mgs.TaxID(i)) for i in range(1, 3)])
        sample_counts = {
            mgs.TaxID(0): {mgs.Sample("a"): 4},
            mgs.TaxID(1): {mgs.Sample("a"): 2, mgs.Sample("b"): 3},
        }
        expected = Counter({mgs.Sample("a"): 6, mgs.Sample("b"): 3})
        self.assertEqual(mgs.count_reads(taxtree, sample_counts), expected)


class TestMGSData(unittest.TestCase):
    # TODO: Fake this?
    mgs_data = MGSData.from_repo()
    bioproject = mgs.BioProject("PRJNA729801")  # Rothman
    sample = mgs.Sample("SRR14530726")  # Random Rothman sample
    taxid = pathogens.pathogens["norovirus"].pathogen_chars.taxid

    def test_from_repo(self):
        self.assertIsInstance(MGSData.from_repo(), MGSData)

    def test_sample_attributes(self):
        samples = self.mgs_data.sample_attributes(self.bioproject)
        self.assertIn(self.sample, samples)
        self.assertIsInstance(samples[self.sample], mgs.SampleAttributes)

    def test_total_reads(self):
        reads = self.mgs_data.total_reads(self.bioproject)
        self.assertIn(self.sample, reads)
        self.assertIsInstance(reads[self.sample], int)

    def test_viral_reads(self):
        reads = self.mgs_data.viral_reads(self.bioproject, self.taxid)
        self.assertIn(self.sample, reads)
        self.assertIsInstance(reads[self.sample], int)


class TestTree(unittest.TestCase):
    leaf = Tree(0)
    node = Tree(0, [Tree(x) for x in range(1, 3)])

    def test_iter(self):
        for i, t in zip(range(1), self.leaf):
            self.assertEqual(i, t.data)
        for i, t in zip(range(3), self.node):
            self.assertEqual(i, t.data)

    def test_get_item(self):
        self.assertEqual(self.leaf, self.leaf[0])
        self.assertIsNone(self.leaf[1])
        self.assertEqual(self.node, self.node[0])
        for i in range(3):
            self.assertIsNotNone(self.node[i])
        self.assertIsNone(self.node[3])

    def test_contains(self):
        self.assertIn(0, self.leaf)
        for i in range(3):
            self.assertIn(i, self.node)

    def test_to_list(self):
        self.assertEqual(self.leaf.to_list(), [0])
        self.assertEqual(self.node.to_list(), [0, [1], [2]])

    def test_parse_inverse(self):
        self.assertEqual(self.node, Tree.tree_from_list(self.node.to_list()))
        self.assertEqual(self.leaf, Tree.tree_from_list(self.leaf.to_list()))

    def test_map(self):
        f = lambda x: x + 1
        self.assertEqual(self.leaf.map(f), Tree(1))
        self.assertEqual(
            self.node.map(f), Tree(f(0), [Tree(f(x)) for x in range(1, 3)])
        )

    def test_map_id(self):
        self.assertEqual(self.node, self.node.map(lambda x: x))
        self.assertEqual(self.leaf, self.leaf.map(lambda x: x))

    def test_map_composition(self):
        f = lambda x: x + 1
        g = lambda x: 2 * x
        self.assertEqual(
            self.node.map(f).map(g), self.node.map(lambda x: g(f(x)))
        )
        self.assertEqual(
            self.leaf.map(f).map(g), self.leaf.map(lambda x: g(f(x)))
        )


if __name__ == "__main__":
    unittest.main()
