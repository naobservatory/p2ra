import dataclasses

from pathogen_properties import *
from populations import us_population

background = """HIV is a sexually-transmitted retrovirus which gradually
weakens the immune system."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(11676),
    selection=SelectionRound.ROUND_1,
)

# Controlled HIV has basically no sheddding, and we're really only interested
# in uncontrolled HIV.  The CDC estimates that 57% of Americans with HIV were
# virally suppressed in 2019.
us_unsuppressed_fraction_2019 = Scalar(
    scalar=1 - 0.57,
    country="United States",
    date="2019",
    source="https://www.cdc.gov/hiv/library/reports/hiv-surveillance/vol-26-no-2/content/national-profile.html#:~:text=57%25%20were%20virally%20suppressed",
)

us_infected_2019 = PrevalenceAbsolute(
    infections=1.2e6,
    country="United States",
    date="2019",
    active=Active.LATENT,
    source="https://www.cdc.gov/hiv/library/reports/hiv-surveillance/vol-26-no-2/content/national-profile.html#:~:text=Among%20the%20estimated-,1.2%20million%20people,-living%20with%20HIV",
)

us_population_2019 = Population(
    people=328_231_337,
    country="United States",
    date="2019-01-01",
    source="https://www.census.gov/newsroom/press-releases/2019/new-years-population.html",
)

la_unsuppressed_fraction_2020 = Scalar(
    scalar=1 - 0.6,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2020",
    source="https://web.archive.org/web/20201202004910/https://www.lacounty.hiv/",
)

la_infected_2020 = PrevalenceAbsolute(
    infections=57_700,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2020",
    active=Active.LATENT,
    source="https://web.archive.org/web/20201202004910/https://www.lacounty.hiv/",
)

denmark_infected_undiagnosed_2022 = PrevalenceAbsolute(
    infections=600,
    country="Denmark",
    date="2022",
    active=Active.LATENT,
    source="https://en.ssi.dk/surveillance-and-preparedness/surveillance-in-denmark/annual-reports-on-disease-incidence/hiv-2020#:~:text=the%20dark%20figure)%3A-,600,-Hereof%20MSM%3A%20300",
    # Matches number from UNAIDS, with a an absolute prevalence of 6_900 and
    # ART coverage of 89% (https://www.unaids.org/en/regionscountries/
    # countries/denmark)
)

hiv_share_copenhagen = Scalar(
    scalar=(
        0.39  # 2017, no raw numbers given
        + ((38 + 16) / (84 + 67))  # 2018
        + ((37 + 8) / (76 + 67))  # 2019
        + ((19 + 9) / (61 + 47))  # 2020
    )
    / 4,
    # Average of 2017- 2020 share of newly diagnosed HIV cases from Copenhagen City, both from HTX and MSM transmission.
    # 2017 data: https://en.ssi.dk/surveillance-and-preparedness/surveillance-in-denmark/annual-reports-on-disease-incidence/hiv-2017
    # 2018 data: https://en.ssi.dk/surveillance-and-preparedness/surveillance-in-denmark/annual-reports-on-disease-incidence/hiv-2018
    # Table 2 in source below shows 2019 and 2020 data.
    source="https://en.ssi.dk/surveillance-and-preparedness/surveillance-in-denmark/annual-reports-on-disease-incidence/hiv-2020#:~:text=rest%20of%20Denmark%2C-,Table%202.,-CD4%20counts%20at",
    date="2020",
)

copenhagen_population_2022 = Population(
    people=652_221,
    country="Denmark",
    date="2022",
    source="https://kk.statistikbank.dk/statbank5a/SelectVarVal/Define.asp?MainTable=KKBEF1"
    # Search variables: Copenhagen in total, Gender in general, Age in total,
    # 2022 Q4
)


def estimate_prevalences() -> list[Prevalence]:
    us_2019 = (
        us_infected_2019.to_rate(us_population_2019)
        * us_unsuppressed_fraction_2019
    )

    # HIV prevalence should be close to constant, so it's fine to figure that
    # it's about the same in 2020 and 2021 as it was in 2019.
    us_2020 = dataclasses.replace(us_2019, date_source=Variable(date="2020"))
    us_2021 = dataclasses.replace(us_2019, date_source=Variable(date="2021"))
    la_2020 = (
        la_infected_2020.to_rate(
            us_population(
                state="California", county="Los Angeles County", year=2020
            )
        )
        * la_unsuppressed_fraction_2020
    )
    # Extrapolating undiagnosed HIV rate in Copenhagen backwards from 2022 to 2015-2018. We assume that HIV has mostly stayed constant.
    copenhagen_2022 = (
        denmark_infected_undiagnosed_2022 * hiv_share_copenhagen
    ).to_rate(copenhagen_population_2022)
    copenhagen_2018 = dataclasses.replace(
        copenhagen_2022, date_source=Variable(date="2018")
    )

    copenhagen_2017 = dataclasses.replace(
        copenhagen_2018, date_source=Variable(date="2017")
    )
    copenhagen_2016 = dataclasses.replace(
        copenhagen_2018, date_source=Variable(date="2016")
    )
    copenhagen_2015 = dataclasses.replace(
        copenhagen_2018, date_source=Variable(date="2015")
    )

    return [
        us_2019,
        us_2020,
        us_2021,
        la_2020,
        copenhagen_2022,
        copenhagen_2018,
        copenhagen_2017,
        copenhagen_2016,
        copenhagen_2015,
    ]


def estimate_incidences() -> list[IncidenceRate]:
    return []
