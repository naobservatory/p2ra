import csv
import datetime
from collections import Counter

from pathogen_properties import *

background = """SARS-CoV-2 is an airborne coronavirus, responsible for the
2019- pandemic"""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(2697049),
)

shedding_duration = SheddingDuration(
    days=18,
    # "In a review of 28 studies, the pooled median duration of viral RNA
    # detection in respiratory specimens was 18 days following the onset of
    # symptoms."
    source="https://www.uptodate.com/contents/covid-19-epidemiology-virology-and-prevention/print",
)

underreporting = Scalar(
    scalar=4.0,
    confidence_interval=(3.4, 4.7),
    coverage_probability=0.95,
    country="United States",
    start_date="2020-02",
    end_date="2021-09",
    source="https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/burden.html#:~:text=1%20in%204.0%20(95%25%20UI*%203.4%20%E2%80%93%204.7)%20COVID%E2%80%9319%20infections%20were%20reported.",
)

county_populations = {
    ("San Diego", "California"): Population(
        people=3_298_635,
        date="2020-04-01",
        country="United States",
        state="California",
        county="San Diego",
        tag="San Diego 2020",
        source="https://www.census.gov/quickfacts/fact/table/sandiegocountycalifornia/PST045221",
    ),
    ("Los Angeles", "California"): Population(
        people=10_014_042,
        date="2020-04-01",
        country="United States",
        state="California",
        county="Los Angeles",
        tag="Los Angeles 2020",
        source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia/PST045222",
    ),
    ("Orange", "California"): Population(
        people=3_186_989,
        date="2020-04-01",
        country="United States",
        state="California",
        county="Orange",
        tag="Orange 2020",
        source="https://www.census.gov/quickfacts/orangecountycalifornia",
    ),
    ("Alameda", "California"): Population(
        people=1_682_353,
        date="2020-04-01",
        country="United States",
        state="California",
        county="Alameda",
        tag="Alameda 2020",
        source="https://www.census.gov/quickfacts/alamedacountycalifornia",
    ),
    ("Marin", "California"): Population(
        people=262_318,
        date="2020-04-01",
        country="United States",
        state="California",
        county="Marin",
        tag="Marin 2020",
        source="https://www.census.gov/quickfacts/marincountycalifornia",
    ),
    ("San Francisco", "California"): Population(
        people=873_959,
        date="2020-04-01",
        country="United States",
        state="California",
        county="San Francisco",
        tag="San Francisco 2020",
        source="https://www.census.gov/quickfacts/sanfranciscocountycalifornia",
    ),
}

ohio_population = Population(
    people=11_799_374,
    date="2020-04-01",
    country="United States",
    state="Ohio",
    tag="Ohio 2020",
    source="https://www.census.gov/quickfacts/OH",
)


def estimate_prevalences():
    estimates = []

    ohio_totals = Counter()  # day -> total new cases, 7d moving average

    # From the COVID-19 Data Repository by the Center for Systems Science and
    # Engineering (CSSE) at Johns Hopkins University
    #
    # Downloaded 2023-05-02 from https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv
    with open(
        prevalence_data_filename("time_series_covid19_confirmed_US.csv")
    ) as inf:
        for row in csv.reader(inf):
            county = row[5]
            state = row[6]

            if state != "Ohio" and (county, state) not in county_populations:
                continue

            # In the tsv file, cumulative case counts start at column 11 with
            # counts for 2020-01-22.
            case_counts = [int(x) for x in row[11:]]
            day = datetime.date.fromisoformat("2020-01-22")

            # For computing a 7-day centered moving average.  We want a moving
            # average because case reporting is not uniform over the week.
            latest = [0] * 7

            for prev_case_count, case_count in zip(
                case_counts, case_counts[1:]
            ):
                # increment day at the beginning because zip means case_count
                # never takes in the initial value.
                day = day + datetime.timedelta(days=1)

                # case counts are cumulative, but we want daily cases
                delta = case_count - prev_case_count
                latest.pop(0)
                latest.append(delta)

                # centered moving average
                # https://www.jefftk.com/p/careful-with-trailing-averages
                date = str(day - datetime.timedelta(days=3))
                annual_infections = sum(latest) * 52
                if state == "Ohio":
                    ohio_totals[date] += annual_infections
                else:
                    cases = IncidenceAbsolute(
                        annual_infections=annual_infections,
                        country="United States",
                        state=state,
                        county=county,
                        date=date,
                        tag="%s 2020" % county,
                    )
                    estimates.append(
                        (
                            cases.to_rate(
                                county_populations[county, state]
                            ).to_prevalence(shedding_duration)
                            * underreporting
                        ).target(
                            country="United States",
                            state=state,
                            county=county,
                            date=date,
                        )
                    )

    for date, annual_infections in ohio_totals.items():
        cases = IncidenceAbsolute(
            annual_infections=annual_infections,
            country="United States",
            state="Ohio",
            date=date,
            tag="Ohio 2020",
        )
        # TODO: we can probably get a better undereporting figure for the
        # omicron surge.
        estimates.append(
            (
                cases.to_rate(ohio_population).to_prevalence(shedding_duration)
                * underreporting
            ).target(
                country="United States",
                state="Ohio",
                date=date,
            )
        )

    return estimates
