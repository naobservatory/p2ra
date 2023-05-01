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


us_national_foodborne_cases_2006 = IncidenceAbsolute(
    annual_infections=5_461_731,
    confidence_interval=(3_227_078, 8_309_480),
    coverage_probability=0.9,  # credible interval
    country="United States",
    tag="us-2006",
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
    tag="us-2006",
    # "all estimates were based on the US population in 2006 (299 million
    # persons)"
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3375761/#:~:text=population%20in%202006%20(-,299%20million%20persons,-).%20Estimates%20were%20derived",
)

shedding_duration = SheddingDuration(
    days=2,
    confidence_interval=(1, 3),
    # "Norovirus infection symptoms usually last 1 to 3 days"
    source="https://www.mayoclinic.org/diseases-conditions/norovirus/symptoms-causes/syc-20355296#:~:text=Norovirus%20infection%20symptoms%20usually%20last%201%20to%203%20days",
)


def days_in_month(date):
    _, last_day = calendar.monthrange(date.year, date.month)
    return last_day


def estimate_prevalences():
    us_outbreaks = Counter()  # date -> count
    ca_outbreaks = Counter()  # date -> count

    prevalences = []

    # Downloaded on 2023-04-28 from https://wwwn.cdc.gov/norsdashboard/
    # Data runs through the end of 2021.
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
            # It distinguishes between GI and GII Norovirus.  I'm currently
            # discarding this, but it could potentially be useful?
            genotype = row[cols.index("Serotype or Genotype")]
            if "Norovirus" not in etiology:
                # It's the National Outbreak Reporting System, not the
                # Norovirus Outbreak Reporting System.
                #
                # The non-Norovirus ones are almost all bacteria or parasites,
                # though, not much useful to us.
                continue

            date = year, month

            us_outbreaks[date] += 1
            if state == "California":
                ca_outbreaks[date] += 1

    normal_year_national_prevalence = (
        us_national_foodborne_cases_2006.to_rate(us_population_2006)
        .to_prevalence(shedding_duration)
        .scale(us_total_relative_to_foodborne_2006)
    )

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
                        date=target_date,
                        source="https://wwwn.cdc.gov/norsdashboard/",
                    )
                ).target(country="United States", date=target_date)
            )

            # The CA-specific ones are very noisy; consider dropping them?
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
