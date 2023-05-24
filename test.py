#!/usr/bin/env python3

import datetime
import unittest
from collections import Counter

import mgs
import pathogens
import populations
import stats
from mgs import MGSData, SampleAttributes
from pathogen_properties import *
from tree import Tree


class TestPathogens(unittest.TestCase):
    def test_hsv1_imported(self):
        self.assertIn("hsv_1", pathogens.pathogens)

    def test_summarize_location(self):
        us_2019, la_2020 = pathogens.pathogens["hiv"].estimate_prevalences()
        self.assertEqual(us_2019.summarize_location(), "United States")
        self.assertEqual(
            la_2020.summarize_location(),
            "Los Angeles County, California, United States",
        )

    def test_dates(self):
        us_2019, la_2020 = pathogens.pathogens["hiv"].estimate_prevalences()
        self.assertEqual(us_2019.parsed_start, datetime.date(2019, 1, 1))
        self.assertEqual(us_2019.parsed_end, datetime.date(2019, 12, 31))

        self.assertEqual(la_2020.parsed_start, datetime.date(2020, 1, 1))
        self.assertEqual(la_2020.parsed_end, datetime.date(2020, 12, 31))

    def test_properties_exist(self):
        for pathogen_name, pathogen in pathogens.pathogens.items():
            with self.subTest(pathogen=pathogen_name):
                self.assertIsInstance(pathogen.background, str)

                self.assertIsInstance(pathogen.pathogen_chars, PathogenChars)

                saw_estimate = False
                for estimate in pathogen.estimate_prevalences():
                    self.assertIsInstance(estimate, Prevalence)
                    saw_estimate = True
                for estimate in pathogen.estimate_incidences():
                    self.assertIsInstance(estimate, IncidenceRate)
                    saw_estimate = True
                self.assertTrue(saw_estimate)

    def test_dates_set(self):
        for pathogen_name, pathogen in pathogens.pathogens.items():
            with self.subTest(pathogen=pathogen_name):
                for estimate in pathogen.estimate_prevalences():
                    estimate.get_dates()


class TestMMWRWeek(unittest.TestCase):
    def test_mmwr_week(self):
        self.assertEqual(
            pathogens.pathogens["influenza"].parse_mmwr_week(2020, 1),
            # Year starts on a Wednesday, so week 1 starts in 2019.
            datetime.date(2019, 12, 29),
        )
        self.assertEqual(
            pathogens.pathogens["influenza"].parse_mmwr_week(2019, 1),
            # Year starts on a Tuesday, so week 1 starts in 2018.
            datetime.date(2018, 12, 30),
        )
        self.assertEqual(
            pathogens.pathogens["influenza"].parse_mmwr_week(2018, 1),
            # Year starts on a Monday, so week 1 starts in 2017.
            datetime.date(2017, 12, 31),
        )
        self.assertEqual(
            pathogens.pathogens["influenza"].parse_mmwr_week(2017, 1),
            # Year starts on a Sunday, so week 1 starts in 2017.
            datetime.date(2017, 1, 1),
        )
        self.assertEqual(
            pathogens.pathogens["influenza"].parse_mmwr_week(2016, 1),
            # Year starts on a Friday, so week 1 starts in 2016.
            datetime.date(2016, 1, 3),
        )
        self.assertEqual(
            pathogens.pathogens["influenza"].parse_mmwr_week(2015, 1),
            # Year starts on a Thursday, so week 1 starts in 2016.
            datetime.date(2015, 1, 4),
        )


class TestVaribles(unittest.TestCase):
    def test_date_parsing(self):
        v = Variable(date="2019")
        self.assertEqual(
            v.get_dates(),
            (datetime.date(2019, 1, 1), datetime.date(2019, 12, 31)),
        )

        v = Variable(date="2019-02")
        self.assertEqual(
            v.get_dates(),
            (datetime.date(2019, 2, 1), datetime.date(2019, 2, 28)),
        )

        v = Variable(date="2020-02")
        self.assertEqual(
            v.get_dates(),
            (datetime.date(2020, 2, 1), datetime.date(2020, 2, 29)),
        )

        v = Variable(date="2020-02-01")
        self.assertEqual(
            v.get_dates(),
            (datetime.date(2020, 2, 1), datetime.date(2020, 2, 1)),
        )

        v = Variable(start_date="2020-01", end_date="2020-02")
        self.assertEqual(
            v.get_dates(),
            (datetime.date(2020, 1, 1), datetime.date(2020, 2, 29)),
        )

        v = Variable(start_date="2020-01-07", end_date="2020-02-06")
        self.assertEqual(
            v.get_dates(),
            (datetime.date(2020, 1, 7), datetime.date(2020, 2, 6)),
        )

        v1 = Variable(date="2019")
        v2 = Variable(date="2020", date_source=v1)
        self.assertEqual(
            v2.get_dates(),
            (datetime.date(2019, 1, 1), datetime.date(2019, 12, 31)),
        )

        with self.assertRaises(Exception):
            Variable(start_date="2020-01-07")

        with self.assertRaises(Exception):
            Variable(end_date="2020-01-07")

        with self.assertRaises(Exception):
            Variable(start_date="2020-01-07", date="2020")

        with self.assertRaises(Exception):
            Variable(end_date="2020-01-07", date="2020")

        with self.assertRaises(Exception):
            Variable(start_date="2020-01-07", end_date="2020-01-06")

        with self.assertRaises(Exception):
            Variable(date="2020-1")

        with self.assertRaises(Exception):
            Variable(date="2020/1/1")

        with self.assertRaises(Exception):
            Variable(date="2020/01/01")

        v = Variable(date="2020")
        with self.assertRaises(AssertionError):
            v.get_date()  # asserts start==end

        v = Variable()
        self.assertIsNone(v.parsed_start)
        self.assertIsNone(v.parsed_end)
        with self.assertRaises(AssertionError):
            v.get_dates()  # asserts dates are set

    def test_locations(self):
        v1 = Variable(
            country="United States", state="Ohio", county="Franklin County"
        )
        self.assertEqual(
            ("United States", "Ohio", "Franklin County"), v1.get_location()
        )

        v2 = Variable(
            country="United States",
            state="California",
            county="Alameda County",
        )

        # Conflicting locations with no resolution specified.
        with self.assertRaises(ValueError):
            Variable(inputs=[v1, v2])

        v3 = Variable(inputs=[v1, v2], location_source=v1)
        self.assertEqual(
            ("United States", "Ohio", "Franklin County"), v3.get_location()
        )


class TestMGS(unittest.TestCase):
    repo = mgs.GitHubRepo(**mgs.MGS_REPO_DEFAULTS)

    def test_load_bioprojects(self):
        bps = mgs.load_bioprojects(self.repo)
        # Rothman
        self.assertIn(mgs.BioProject("PRJNA729801"), bps)

    def test_load_sample_attributes(self):
        samples = mgs.load_sample_attributes(self.repo)
        # Randomly picked Rothman sample
        s = mgs.Sample("SRR14530726")
        self.assertIn(s, samples)
        attrs = samples[s]
        self.assertEqual(attrs.country, "United States")
        self.assertEqual(attrs.state, "California")
        self.assertEqual(attrs.county, "San Diego County")
        self.assertEqual(attrs.date, datetime.date(2020, 8, 27))
        self.assertEqual(attrs.enrichment, mgs.Enrichment.VIRAL)

    def test_load_sample_counts(self):
        sample_counts = mgs.load_sample_counts(self.repo)
        for p in ["sars_cov_2", "hiv"]:
            with self.subTest(pathogen=p):
                for taxid in pathogens.pathogens[p].pathogen_chars.taxids:
                    self.assertIn(taxid, sample_counts)

    def test_load_tax_tree(self):
        tree = mgs.load_tax_tree(self.repo)
        for pathogen_name, pathogen in pathogens.pathogens.items():
            with self.subTest(pathogen=pathogen_name):
                for taxid in pathogen.pathogen_chars.taxids:
                    self.assertIn(taxid, tree)

    def test_count_reads(self):
        taxtree = Tree(mgs.TaxID(0), [Tree(mgs.TaxID(i)) for i in range(1, 3)])
        sample_counts = {
            mgs.TaxID(0): {mgs.Sample("a"): 4},
            mgs.TaxID(1): {mgs.Sample("a"): 2, mgs.Sample("b"): 3},
        }
        expected = Counter({mgs.Sample("a"): 6, mgs.Sample("b"): 3})
        self.assertEqual(mgs.count_reads(taxtree, sample_counts), expected)


class TestMGSData(unittest.TestCase):
    mgs_data = MGSData.from_repo()
    bioproject = mgs.BioProject("PRJNA729801")  # Rothman
    sample = mgs.Sample("SRR14530726")  # Random Rothman sample
    taxids = pathogens.pathogens["norovirus"].pathogen_chars.taxids

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
        reads = self.mgs_data.viral_reads(self.bioproject, self.taxids)
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


class TestPopulations(unittest.TestCase):
    def test_county_state(self):
        self.assertEqual(
            populations.us_population(
                county="Bristol County", state="Rhode Island", year=2020
            ),
            Population(
                people=50_774,
                date="2020-07-01",
                source="https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html",
                country="United States",
                state="Rhode Island",
                county="Bristol County",
            ),
        )

        self.assertEqual(
            populations.us_population(
                county="Bristol County", state="Rhode Island", year=2020
            ).people,
            50_774,
        )
        self.assertEqual(
            populations.us_population(
                county="Bristol County", state="Rhode Island", year=2021
            ).people,
            50_800,
        )
        self.assertEqual(
            populations.us_population(
                county="Bristol County", state="Rhode Island", year=2022
            ).people,
            50_360,
        )

        self.assertEqual(
            populations.us_population(
                county="Southeastern Connecticut Planning Region",
                state="Connecticut",
                year=2022,
            ).people,
            280_403,
        )

    def test_state(self):
        self.assertEqual(
            populations.us_population(state="California", year=2022),
            Population(
                # From https://www.census.gov/quickfacts/CA
                people=39_029_342,
                date="2022-07-01",
                source="https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html",
                country="United States",
                state="California",
            ),
        )

    def test_country(self):
        self.assertEqual(
            populations.us_population(year=2022),
            Population(
                # https://www.census.gov/quickfacts/USA
                people=333_287_557,
                date="2022-07-01",
                source="https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html",
                country="United States",
            ),
        )


class TestStats(unittest.TestCase):
    attrs = SampleAttributes(
        country="United States",
        state="Pennsylvania",
        county="Allegheny County",
        date=datetime.date.fromisoformat("2019-05-14"),
        reads=100,
        location="Loc",
    )

    def test_is_match(self):
        v1 = Variable(country="United States", date="2019")
        self.assertTrue(stats.is_match(self.attrs, v1))
        v2 = Variable(country="United States", date="2019-05-14")
        self.assertTrue(stats.is_match(self.attrs, v2))
        v3 = Variable(country="United States", date="2019-05-15")
        self.assertFalse(stats.is_match(self.attrs, v3))
        v4 = Variable(
            country="United States",
            start_date="2019-05-01",
            end_date="2019-06-02",
        )
        self.assertTrue(stats.is_match(self.attrs, v4))
        v5 = Variable(
            country="United States",
            state="Pennsylvania",
            county="Allegheny County",
            date="2019",
        )
        self.assertTrue(stats.is_match(self.attrs, v5))
        v6 = Variable(
            country="United States",
            state="Pennsylvania",
            county="Beaver County",
            date="2019",
        )
        self.assertFalse(stats.is_match(self.attrs, v6))
        v7 = Variable(
            country="United States",
            state="Ohio",
            county="Lake County",
            date="2019",
        )
        self.assertFalse(stats.is_match(self.attrs, v7))

    def test_lookup_variable(self):
        v1 = Variable(country="United States", date="2019")
        v2 = Variable(country="United States", date="2019-05-14")
        v3 = Variable(country="United States", date="2019-05-15")
        self.assertEqual(stats.lookup_variable(self.attrs, [v1, v3]), v1)
        self.assertEqual(stats.lookup_variable(self.attrs, [v2, v3]), v2)
        self.assertEqual(stats.lookup_variable(self.attrs, [v3, v1]), v1)
        with self.assertRaises(AssertionError):
            stats.lookup_variable(self.attrs, [v1, v2])
        with self.assertRaises(AssertionError):
            stats.lookup_variable(self.attrs, [v3])


if __name__ == "__main__":
    unittest.main()
