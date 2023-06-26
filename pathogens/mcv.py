import dataclasses

from pathogen_properties import *

background = """Merkel cell polyomavirus is a common virus that is
suspected to occasionally cause skin cancer."""
pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(493803),
)

seroprevalence_women_us_2009 = Prevalence(
    infections_per_100k=59.4 * 100_000,
    number_of_participants=451,
    country="United States",
    date="2009",
    active=Active.LATENT,
    source="https://academic.oup.com/jnci/article/101/21/1510/963538#:~:text=A%20cut%20point%20of%2015%E2%80%89000%20median%20fluorescent%20intensity%20units%20was%20chosen%20as%20described%20above%20and%20used%20to%20identify%20268%20(59.4%25)%20of%20the%20451%20control%20subjects%20as%20seropositive%20for%20MCPyV%20VP1%3B%20this%20percentage%20was%20similar%20to%20that%20among%20women%20in%20control%20group%201%20(ie%2C%2053%25)."
    # Though this seroprevalence number is only for women, it tracks 2011 data
    # from Italy, showing a seroprevalence of around 70% in a control
    # population of 945 individuals. Table 1, plotted in Figure 4:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3187023/#:~:text=MCPyV%2C%20JCPyV%2C%20and%20BKPyV%20age%2Dspecific%20seroprevalencea
)


def estimate_prevalences() -> list[Prevalence]:
    # There is no vaccine against MCV, so seroprevalence should stay
    # approximately constant.
    us_2020 = dataclasses.replace(
        seroprevalence_women_us_2009, date_source=Variable(date="2020")
    )
    us_2021 = dataclasses.replace(
        seroprevalence_women_us_2009, date_source=Variable(date="2021")
    )

    return [
        us_2020,
        us_2021,
    ]


def estimate_incidences():
    return []
