from pathogen_properties import *

background = """JC polyomavirus is a polyomavirus that is ubiquitous in the human population, infecting 70-90% of people worldwide and staying latent in the body thereafter. It is the causative agent of progressive multifocal leukoencephalopathy (PML), a rare but often fatal demyelinating disease of the central nervous system."""

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(10632),
)

ch_2009_seroprevalence = Prevalence(
    infections_per_100k=0.58 * 100_000,
    # "JCV IgG seroprevalence of 58% (231 of 400)"
    number_of_participants=400,
    # Note that this study also provides numbers on urinary shedding:
    # "JCV, [...] was detected in samples from 75 (19%) of the donors"
    # Given that JC Virus stays latent after infection, seroprevalence likely
    # tracks prevalence more closely than urinary shedding.
    country="Switzerland",
    date="2009",
    active=Active.LATENT,
    source="https://academic.oup.com/jid/article/199/6/837/2192120?login=false#90076999:~:text=Prevalence%20of%20BKV%20and%20JCV%20antibodies%2C%20DNAemia%2C%20and%20DNAuriaBKV%20IgG%20seroprevalence%20among%20healthy%20blood%20donors%20was%2082%25%20(328%20of%20400)%2C%20significantly%20higher%20than%20the%20corresponding%20JCV%20IgG%20seroprevalence%20of%2058%25%20(231%20of%20400)",
)


ch_2009_urine_shedding
