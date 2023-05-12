import csv
import datetime
from collections import Counter

from pathogen_properties import *
from populations import us_population

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

target_counties = set(
    [
        "San Diego County, California",
        "Los Angeles County, California",
        "Orange County, California",
        "Alameda County, California",
        "Marin County, California",
        "San Francisco County, California",
        "Franklin County, Ohio",
        "Greene County, Ohio",
        "Lawrence County, Ohio",
        "Licking County, Ohio",
        "Lucas County, Ohio",
        "Montogmery County, Ohio",
        "Sandusky County, Ohio",
        "Summit County, Ohio",
        "Trumbull County, Ohio",
    ]
)


def estimate_prevalences():
    estimates = []

    # From the COVID-19 Data Repository by the Center for Systems Science and
    # Engineering (CSSE) at Johns Hopkins University
    #
    # Downloaded 2023-05-02 from https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv
    incidences = []

    with open(
        prevalence_data_filename("time_series_covid19_confirmed_US.csv")
    ) as inf:
        for row in csv.reader(inf):
            truncated_county = row[5]
            state = row[6]

            # The time series data names counties like "Franklin", but the full
            # name the census uses is "Franklin County".  Use the full names
            # for consistency.
            county = "%s County" % truncated_county
            county_state = "%s, %s" % (county, state)

            if county_state not in target_counties:
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
                date = day - datetime.timedelta(days=3)
                if date.year > 2022:
                    continue

                # Right now we use the same underreporting figure for both
                # Spring/Fall 2020 and Winter 2021-2022.
                #
                # TODO: we can probably get a better undereporting figure for
                # the omicron surge and this is likely too small.  The CDC 4x
                # figure is not intended to cover this time period, this was
                # after rapid tests were starting to be available, and omicron
                # was relatively mild.

                annual_infections = sum(latest) * 52

                cases = IncidenceAbsolute(
                    annual_infections=annual_infections,
                    country="United States",
                    state=state,
                    county=county,
                    date=date.isoformat(),
                )
                incidences.append(
                    cases.to_rate(
                        us_population(
                            county=county, state=state, year=date.year
                        )
                    ).target(
                    * underreporting
                )

    return shedding_duration.prevalences_from_incidences(incidences)
