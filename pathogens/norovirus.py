from collections import Counter, defaultdict

import numpy as np

from pathogen_properties import *

background = """Norovirus in a GI infection, mostly spread through personal
contact."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(142786),
)


normal_year_national_cases = IncidenceAbsolute(
    annual_infections=20e6,
    confidence_interval=(19e6, 21e6),
    country="United States",
    tag="us",
    source="https://www.cdc.gov/norovirus/trends-outbreaks/burden-US.html",
)

us_population = Population(
    people=333_287_557,
    country="United States",
    date="2022-07-01",
    tag="us",
    source="https://www.census.gov/quickfacts/fact/table/US/PST045221",
)

shedding_duration = SheddingDuration(
    days=2,
    confidence_interval=(1, 3),
    source="https://www.mayoclinic.org/diseases-conditions/norovirus/symptoms-causes/syc-20355296",
)


def days_in_month(date):
    _, last_day = calendar.monthrange(date.year, date.month)
    return last_day


def estimate_prevalences():
    us_outbreaks = Counter()  # date -> count
    ca_outbreaks = Counter()  # date -> count

    prevalences = []

    # Downloaded from https://wwwn.cdc.gov/norsdashboard/ 2023-04-28.
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
            etiology = row[cols.index("Etiology")]
            genotype = row[cols.index("Serotype or Genotype")]
            if "Norovirus" not in etiology:
                continue

            date = year, month

            us_outbreaks[date] += 1
            if state == "California":
                ca_outbreaks[date] += 1

    normal_year_national_prevalence = normal_year_national_cases.to_rate(
        us_population
    ).to_prevalence(shedding_duration)

    total_us_outbreaks = 0
    total_ca_outbreaks = 0
    days_considered = 0
    for year in range(2012, 2020):
        for month in range(1, 13):
            total_us_outbreaks += us_outbreaks[year, month]
            total_ca_outbreaks += ca_outbreaks[year, month]
            days_considered += days_in_month(datetime.date(year, month, 1))
    normal_us_average_daily_outbreaks = total_us_outbreaks / days_considered
    normal_ca_average_daily_outbreaks = total_ca_outbreaks / days_considered

    for target_year in range(2012, 2022):
        for target_month in range(1, 13):
            target_date = "%s-%s" % (target_year, str(target_month).zfill(2))

            target_us_daily_outbreaks = us_outbreaks[
                target_year, target_month
            ] / days_in_month(datetime.date(target_year, target_month, 1))
            target_ca_daily_outbreaks = ca_outbreaks[
                target_year, target_month
            ] / days_in_month(datetime.date(target_year, target_month, 1))

            prevalences.append(
                normal_year_national_prevalence.scale(
                    Scalar(
                        scalar=target_us_daily_outbreaks
                        / normal_us_average_daily_outbreaks,
                        country="United States",
                        state="California",
                        date=target_date,
                        source="https://wwwn.cdc.gov/norsdashboard/",
                    )
                ).target(country="United States", date=target_date)
            )

            prevalences.append(
                normal_year_national_prevalence.scale(
                    Scalar(
                        scalar=target_ca_daily_outbreaks
                        / normal_ca_average_daily_outbreaks,
                        country="United States",
                        state="California",
                        date=target_date,
                        source="https://wwwn.cdc.gov/norsdashboard/",
                    )
                ).target(
                    country="United States",
                    state="California",
                    date=target_date,
                )
            )
    return prevalences
