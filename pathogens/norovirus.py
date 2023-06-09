import dataclasses
from collections import Counter, defaultdict

import numpy as np

from pathogen_properties import *

background = """Norovirus is a GI infection, mostly spread through personal
contact."""

NOROVIRUS_GROUP_I = TaxID(122928)
NOROVIRUS_GROUP_II = TaxID(122929)

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxids=frozenset((NOROVIRUS_GROUP_I, NOROVIRUS_GROUP_II)),
    names_by_taxid={
        NOROVIRUS_GROUP_I: "Norovirus (GI)",
        NOROVIRUS_GROUP_II: "Norovirus (GII)",
    },
    selection=SelectionRound.ROUND_1,
)

# We're using Scallan 2011
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3375761/) which:
#
#  1. Estimates the annual number of foodborne Norovirus cases
#  2. Estimates the fractionof Norovirus cases that are foodborne
#
# This gives us enough information to estimate the annual number of Norovirus
# cases.
#
# This is the source for the CDC's estimate:
# https://wwwnc.cdc.gov/eid/article/19/8/pdfs/13-0465.pdf cites this paper as
# their source for 21M annual cases in the US.
#
# Note that this is "2006" as in "relative to the 2006 population" not as in
# "number of cases in 2006".  They're using several years worth of data
# ("mostly from 2000–2008") to make their estimates.

us_national_foodborne_cases_2006 = IncidenceAbsolute(
    annual_infections=5_461_731,
    confidence_interval=(3_227_078, 8_309_480),
    coverage_probability=0.9,  # credible interval
    country="United States",
    date="2006",
    # "Domestically acquired foodborne, mean (90% credible interval)
    # ... 5,461,731 (3,227,078–8,309,480)"
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3375761/#:~:text=5%2C461%2C731%20(3%2C227%2C078%E2%80%938%2C309%2C480)",
)

us_total_relative_to_foodborne_2006 = Scalar(
    scalar=1 / 0.26,
    country="United States",
    date="2006",
    # "Foodborne % ... 26"
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3375761/#:~:text=%3C1-,26,-5%2C461%2C731%20(3%2C227%2C078%E2%80%938%2C309%2C480",
)

us_population_2006 = Population(
    people=299_000_000,
    country="United States",
    date="2006",
    # "all estimates were based on the US population in 2006 (299 million
    # persons)"
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3375761/#:~:text=population%20in%202006%20(-,299%20million%20persons,-).%20Estimates%20were%20derived",
)

monthwise_count = dict[tuple[int, int], float]  # [year, month] -> float


def to_daily_counts(outbreaks_per_month: monthwise_count) -> monthwise_count:
    outbreaks_per_day: monthwise_count = {}
    for year, month in outbreaks_per_month:
        outbreaks_per_day[year, month] = outbreaks_per_month[
            year, month
        ] / days_in_month(year, month)
    return outbreaks_per_day


def load_nors_outbreaks() -> (
    tuple[monthwise_count, monthwise_count, monthwise_count]
):
    us_outbreaks: monthwise_count = defaultdict(float)  # date -> count
    us_outbreaks_I: monthwise_count = defaultdict(float)  # date -> count
    us_outbreaks_II: monthwise_count = defaultdict(float)  # date -> count

    # seen_I, seen_II -> count
    totals: dict[tuple[bool, bool], int] = Counter()
    total_seen_other = 0

    # Downloaded on 2023-04-28 from https://wwwn.cdc.gov/norsdashboard/
    # Click "Download all NORS Dashboard data (Excel)."
    # Exported from Google Sheets as CSV.
    #
    # Data entries start in 1971 and run through the end of 2021.
    # We're only using data from HISTORY_START (2012) onward.
    with open(prevalence_data_filename("cdc-nors-outbreak-data.tsv")) as inf:
        cols = None
        for line in inf:
            row = line.strip().split("\t")
            if not cols:
                cols = row
                continue

            year = int(row[cols.index("Year")])
            month = int(row[cols.index("Month")])
            state = row[cols.index("State")]
            etiologies = row[cols.index("Etiology")]
            # It distinguishes between GI and GII Norovirus.  I'm currently
            # discarding this, but it could potentially be useful?
            genotype = row[cols.index("Serotype or Genotype")]
            if "Norovirus" not in etiologies:
                # It's the National Outbreak Reporting System, not the
                # Norovirus Outbreak Reporting System.
                #
                # The non-Norovirus ones are almost all bacteria or parasites,
                # though, not much useful to us.
                continue

            date = year, month

            seen_I = False
            seen_II = False
            seen_other = False
            for etiology in etiologies.split("; "):
                if etiology.endswith("Norovirus Genogroup I"):
                    seen_I = True
                elif etiology.endswith("Norovirus Genogroup II"):
                    seen_II = True
                elif "Genogroup" in etiology:
                    seen_other = True

            us_outbreaks[date] += 1
            if seen_I and not seen_II:
                us_outbreaks_I[date] += 1
            elif seen_II and not seen_I:
                us_outbreaks_II[date] += 1

            if seen_other:
                total_seen_other += 1

            # We don't care about the I-vs-II labeling in old data, so ignore
            # dates before HISTORY_START.
            if year >= HISTORY_START:
                totals[seen_I, seen_II] += 1

    total_classified = (
        totals[True, False] + totals[False, True] + totals[True, True]
    )

    seen_both_fraction = totals[True, True] / total_classified

    # As of 2023-05-03 this was 1.05%, low enough to ignore.  If this were
    # higher we'd need to estimate prevalences that didn't add to the total
    # prevalence.
    assert seen_both_fraction < 0.011

    # As of 2023-05-03 this was 0.15%, low enough to ignore.  If this were
    # non-trivial we might want to try assigning reads to other genogroups.
    assert total_seen_other / total_classified < 0.0015

    return us_outbreaks, us_outbreaks_I, us_outbreaks_II


def determine_average_daily_outbreaks(us_outbreaks: monthwise_count) -> float:
    total_us_outbreaks = 0.0
    days_considered = 0
    for year in range(HISTORY_START, COVID_START):
        for month in range(1, 13):
            total_us_outbreaks += us_outbreaks[year, month]
            days_considered += days_in_month(year, month)
    return total_us_outbreaks / days_considered


# When estimating the historical pattern, use 2012 through 2019.  This is:
#  * Recent enough to have good data
#  * Long enough to reduce noise
#  * Pre-covid
HISTORY_START = 2012
COVID_START = 2020
DATA_END = 2022


def estimate_incidences() -> list[IncidenceRate]:
    incidences = []

    us_outbreaks, us_outbreaks_I, us_outbreaks_II = load_nors_outbreaks()
    pre_covid_us_average_daily_outbreaks = determine_average_daily_outbreaks(
        us_outbreaks
    )

    pre_covid_national_incidence = (
        us_national_foodborne_cases_2006.to_rate(us_population_2006)
        * us_total_relative_to_foodborne_2006
    )

    us_daily_outbreaks = to_daily_counts(us_outbreaks)

    for year in range(HISTORY_START, DATA_END):
        us_total_I = 0.0
        us_total_II = 0.0
        for month in range(1, 13):
            us_total_I += us_outbreaks_I[year, month]
            us_total_II += us_outbreaks_II[year, month]
        assert us_total_I
        assert us_total_II

        for month in range(1, 13):
            adjustment = Scalar(
                scalar=us_daily_outbreaks[year, month]
                / pre_covid_us_average_daily_outbreaks,
                country="United States",
                date=f"{year}-{month:02d}",
                source="https://wwwn.cdc.gov/norsdashboard/",
            )
            adjusted_national_incidence = dataclasses.replace(
                pre_covid_national_incidence * adjustment,
                date_source=adjustment,
            )

            # Assume that all Norovirus infections are either Group I or II,
            # which is very close (see assertion above).  Also assume the
            # outbreaks for which we have subtype info are representative of
            # all infections.
            us_I = us_outbreaks_I[year, month]
            us_II = us_outbreaks_II[year, month]
            if us_I and us_II:
                group_I_fraction = us_I / (us_I + us_II)
                group_II_fraction = us_II / (us_I + us_II)
            else:
                # Not enough records for this month to compute a ratio between
                # I and II.  Fall back to the annual ratio.
                group_I_fraction = us_total_I / (us_total_I + us_total_II)
                group_II_fraction = us_total_II / (us_total_I + us_total_II)

            incidences.append(
                dataclasses.replace(
                    adjusted_national_incidence
                    * Scalar(scalar=group_I_fraction),
                    taxid=NOROVIRUS_GROUP_I,
                )
            )
            incidences.append(
                dataclasses.replace(
                    adjusted_national_incidence
                    * Scalar(scalar=group_II_fraction),
                    taxid=NOROVIRUS_GROUP_II,
                )
            )

    return incidences


def estimate_prevalences() -> list[Prevalence]:
    return []
