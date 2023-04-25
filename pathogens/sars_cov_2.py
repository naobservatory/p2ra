import datetime

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
    country="United States",
    start_date="2020-02",
    end_date="2021-09",
    source="https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/burden.html#:~:text=1%20in%204.0%20(95%25%20UI*%203.4%20%E2%80%93%204.7)%20COVID%E2%80%9319%20infections%20were%20reported.",
)

county_populations = {
    "San Diego": Population(
        people=3_298_635,
        date="2020-04-01",
        country="United States",
        state="California",
        county="San Diego",
        tag="San Diego 2020",
        source="https://www.census.gov/quickfacts/fact/table/sandiegocountycalifornia/PST045221",
    ),
    "Los Angeles": Population(
        people=10_014_042,
        date="2020-04-01",
        country="United States",
        state="California",
        county="Los Angeles",
        tag="Los Angeles 2020",
        source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia/PST045222",
    ),
    "Orange": Population(
        people=3_186_989,
        date="2020-04-01",
        country="United States",
        state="California",
        county="Orange",
        tag="Orange 2020",
        source="https://www.census.gov/quickfacts/orangecountycalifornia",
    ),
}


def estimate_prevalences():
    estimates = []

    with open(
        prevalence_data_filename("southern-ca-cumulative-covid-cases.tsv")
    ) as inf:
        for line in inf:
            bits = line.strip().split("\t")
            county = bits[5]

            # In the tsv file, cumulative case counts start at column 11 with
            # counts for 2020-01-22.
            case_counts = [int(x) for x in bits[11:]]
            day = datetime.date.fromisoformat("2020-01-22")

            # For computing a 7-day centered moving average.  We want a moving
            # average because case reporting is not uniform over the week.
            latest = [0] * 7

            for prev_case_count, case_count in zip(
                case_counts, case_counts[1:]
            ):
                # increment day at the beginning because zip means case_count never
                # takes in the initial value.
                day = day + datetime.timedelta(days=1)

                # case counts are cumulative, but we want daily cases
                delta = case_count - prev_case_count
                latest.pop(0)
                latest.append(delta)

                # centered moving average
                # https://www.jefftk.com/p/careful-with-trailing-averages
                date = str(day - datetime.timedelta(days=3))
                cases = IncidenceAbsolute(
                    annual_infections=sum(latest) * 52,
                    country="United States",
                    state="California",
                    county=county,
                    date=date,
                    tag="%s 2020" % county,
                )
                estimates.append(
                    cases.to_rate(county_populations[county])
                    .to_prevalence(shedding_duration)
                    .scale(underreporting)
                    .target(
                        country="United States",
                        state="California",
                        county=county,
                        date=date,
                    )
                )

        return estimates
