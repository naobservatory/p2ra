import dataclasses

from pathogen_properties import *

background = """Merkel cell polyomavirus is a common virus that is
suspected to occasionally cause skin cancer."""
pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(493803),
    selection=SelectionRound.ROUND_2,
)


nl_seroprevalence_2018 = Prevalence(
    infections_per_100k=0.819 * 100_000,
    # Found in Table 2, row 5, column 1.
    # Note that we didn't use this study for JCV and BKV, as the underlying
    # methods paper for the test used shows that the two tests used for JCV
    # and BKV show some crossreactivity and are thus likely to overestimate
    # seroprevalence: https://journals.asm.org/doi/full/10.1128/jcm.01566-17#:~:text=Preincubation%20with%20JCPyV,S6B1%20and%20B3
    country="Netherlands",
    number_of_participants=1044,
    # The study population consisted of serum samples from 1050 Dutch blood
    # donors. Donors were included using weighted random selection from Dutch
    # blood donations to obtain groups of equal size in terms of age and sex.
    # Serum samples from eighty blood donation centres were collected over a
    # period of two weeks to ensure an even geographic distribution over the
    # Netherlands.
    date="2018",
    active=Active.LATENT,
    source="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0206273#pone.0206273.ref024:~:text=Table%202.%20Seropositivity%20numbers%20and%20seroprevalence.",
)

seroprevalence_women_us_2009 = Prevalence(
    infections_per_100k=0.594 * 100_000,
    number_of_participants=451,
    country="United States",
    date="2009",
    active=Active.LATENT,
    source="https://academic.oup.com/jnci/article/101/21/1510/963538#:~:text=A%20cut%20point%20of%2015%E2%80%89000%20median%20fluorescent%20intensity%20units%20was%20chosen%20as%20described%20above%20and%20used%20to%20identify%20268%20(59.4%25)%20of%20the%20451%20control%20subjects%20as%20seropositive%20for%20MCPyV%20VP1%3B%20this%20percentage%20was%20similar%20to%20that%20among%20women%20in%20control%20group%201%20(ie%2C%2053%25)."
    # Other studies give numbers that are up to 20 percentage points smaller
    # or larger:
    # "In an early study from the United States (US), MCPyV seroprevalence in
    # adults ranged between 25% (MCPyV strain 350) and 42% (major capsid
    # protein VP1 of MCPyV strain 339) [11]. In more recent studies, higher
    # MCPyV (VP1) seroprevalence rates were found: 82% in blood donors from
    # the Netherlands [17], 69% in adults from Hungary [18], and up to 96% in
    # Italian persons aged 70–79 years"
    # Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9776808/#:~:text=In%20an%20early,70%E2%80%9379%20years
    # Based on this, our number is at least not an outlier, so even though it
    # only measures female seroprevalence, we'll use it.
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
    # Due to a lack of polyomavirus prevalence data for Denmark, we
    # extrapolate Dutch data from 2018 to Denmark, 2015-2018.
    dk_2015 = dataclasses.replace(
        nl_seroprevalence_2018,
        date_source=Variable(date="2015"),
        location_source=Variable(country="Denmark"),
    )
    dk_2016 = dataclasses.replace(
        nl_seroprevalence_2018,
        date_source=Variable(date="2016"),
        location_source=Variable(country="Denmark"),
    )
    dk_2017 = dataclasses.replace(
        nl_seroprevalence_2018,
        date_source=Variable(date="2017"),
        location_source=Variable(country="Denmark"),
    )
    dk_2018 = dataclasses.replace(
        nl_seroprevalence_2018,
        date_source=Variable(date="2018"),
        location_source=Variable(country="Denmark"),
    )
    return [
        us_2020,
        us_2021,
        nl_seroprevalence_2018,
        dk_2015,
        dk_2016,
        dk_2017,
        dk_2018,
    ]


def estimate_incidences():
    return []
