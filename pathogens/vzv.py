from pathogen_properties import *

background = """Varicella-zoster virus (VZV) is a double-stranded DNA virus 
that belongs to the herpesvirus family. It is the virus responsible for 
causing two distinct clinical conditions: varicella (commonly known as 
chickenpox) and herpes zoster (commonly known as shingles).

After a primary infection with VZV, the virus can remain dormant in the 
sensory nerve cells of the body. Later in life, VZV can reactivate and cause 
herpes zoster, or shingles. Herpes zoster is characterized by a painful, 
blistering rash that usually occurs in a localized area along a specific nerve
pathway. Herpes zoster typically occurs in older adults or individuals with 
weakened immune systems, although it can occur in people of any age.

There is a vaccine for VZV. The VZV vaccine became available in the US in 1995 
(https://www.cdc.gov/vaccines/vpd/varicella/public/index.html#:~:text=Chickenpox%20vaccine%20became%20available%20in%20the%20United%20States%20in%201995.). 
This was around 28 years ago."""


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10335),
)

years_since_vaccine = Scalar(
    scalar=2023 - 1995,
    source="https://www.cdc.gov/vaccines/vpd/varicella/public/index.html#:~:text=Chickenpox%20vaccine%20became%20available%20in%20the%20United%20States%20in%201995.",
)

""" The VZV vaccine became available in 1995, 28 years ago. According to the 
CDC, in the prevaccine era, around 40% of people got chickenpox before age 4, 
and around 90% before age 15. From these numbers, I estimate that around  half 
the population got VZV before age 6. So, I'm using the fraction of the US 
population 34 and older as an estimate for the proportion of the population 
that got VZV in the prevaccine era. The median age in the US is 38, so around 
50% of people are 38+. I estimate that around 55% are 34+. 
https://www.cdc.gov/vaccines/pubs/pinkbook/varicella.html#:~:text=In%20the%20prevaccine%20era%2C%20varicella,younger%20than%20age%2015%20years.
https://ourworldindata.org/grapher/median-age"""
population_above_34 = Prevalence(
    infections_per_100k=0.55 * 100000,
    date="2021",
    active=Active.LATENT,
    source="https://ourworldindata.org/grapher/median-age",
)


annual_new_VZV_infections = Scalar(
    scalar=4000000,
    source="https://www.childrenshospital.org/conditions/chickenpox#:~:text=More%20than%2095%20percent%20of%20American%20adults%20have%20had%20chickenpox%20and%20about%204%2C000%2C000%20people%20get%20chickenpox%20every%20year",
)

US_Population = Population(
    people=335000000,
    tag="US-2023",
    source="https://www.census.gov/popclock/",
)

prevalence_from_postvaccine_cases = Prevalence(
    infections_per_100k=annual_new_VZV_infections.scalar
    / US_Population.people
    * 100000
    * years_since_vaccine.scalar,
    active=Active.LATENT,
)

# prevaccine cases estimate + postvaccine cases estimate
absolute_latent_prevalence = Prevalence(
    infections_per_100k=population_above_34.infections_per_100k
    + prevalence_from_postvaccine_cases.infections_per_100k,
    active=Active.LATENT,
)


annual_new_shingles_infections = Scalar(
    scalar=1000000,
    source="https://www.cdc.gov/shingles/hcp/clinical-overview.html#:~:text=An%20estimated%20one%20million%20cases%20of%20herpes%20zoster%20occur%20annually%20in%20the%20United%20States",
)

active_prevalence = PrevalenceAbsolute(
    infections=annual_new_shingles_infections * shingles_shedding_duration
    + annual_new_VZV_infections * chickenpox_shedding_duration,
)


def estimate_prevalences():
    return [absolute_latent_prevalence]
