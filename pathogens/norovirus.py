from pathogen_properties import *

background = """Norovirus in a GI infection, mostly spread through personal
contact."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(142786),
)


cases = IncidenceAbsolute(
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

rothman_period_outbreaks = Number(
    number=13 / 5,  # 13 outbreaks over 5 months
    country="United States",
    source="https://www.cdc.gov/norovirus/reporting/norostat/data-table.html",
    start_date="2020-08-01",
    end_date="2020-12-31",
)
normal_year_outbreaks = Number(
    number=1246 / 12,  # 1246 outbreaks over one year
    country="United States",
    source="https://www.cdc.gov/norovirus/reporting/norostat/data-table.html",
    start_date="2012",
    end_date="2022",
)


def estimate_prevalences():
    return [
        cases.to_rate(us_population)
        .to_prevalence(shedding_duration)
        .scale(rothman_period_outbreaks.per(normal_year_outbreaks))
    ]
