#!/usr/bin/env python3
import sys

sys.path.append("..")

from p2ra_prevalences.pathogen_properties import *
import math

background = """Hepatitis A is a vaccine-preventable liver-infection caused by the hepatitis A virus. In the US,
it's mostly spread by individual contact. There is little seasonal variance in Hepatitis A incidence (https://www.cdc.gov/vaccines/pubs/pinkbook/hepa.html#:~:text=There%20is%20no%20appreciable%20seasonal%20variation%20in%20hepatitis%20A%20incidence.%20In%20the%20prevaccine%20era%2C%20cyclic%20increases%20in%20reported%20acute%20cases%20were%20observed%20every%2010%20to%2015%20years%20and%20were%20characterized%20by%20large%20community%20outbreaks%20of%20disease.).
Viral shedding persists for 1 to 3 weeks."""

pathogen_chars = PathogenChars(
    na_type=NAType.RNA.value,
    enveloped=Enveloped.NON_ENVELOPED.value,
    taxid=208726,
)

prevalence_estimators = {
    "us_incidence_2018": PrevalenceEstimator(
        value: 12,474,
        value_type: yearly_incidence,
        

prevalence_vars = {
    "us_incidence_rate_2018": PrevalenceVariable(
        variable_type=VariableType.MEASUREMENT.value,
        incidence_rate=3.8,
        country="United States",
        source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/index.htm",
    ),
    "us_average_prevalence": PrevalenceVariable(
        variable_type=VariableType.CALCULATED.value,
        prevalence=0.14,
        country="United States",
        source="[please fill in manually]",
    ),
    "washington_state_incidence_rate": PrevalenceVariable(
        variable_type=VariableType.MEASUREMENT.value,
        incidence_rate=[0.4, 0.4, 0.5],
        country="United States",
        state="Washington",
        source="https://www.doh.wa.gov/DataandStatisticalReports/DiseasesandChronicConditions/HepatitisA",
    ),
    "king_county_incidence_rate": PrevalenceVariable(
        variable_type=VariableType.MEASUREMENT.value,
        incidence_rate=[0.5, 0.6],
        country="United States",
        state="Washington",
        county="King",
        source="[please fill in manually]",
    ),
    "king_county_prevalence": PrevalenceVariable(
        variable_type=VariableType.CALCULATED.value,
        prevalence=0.04,
        country="United States",
        state="Washington",
        county="King",
        source="[please fill in manually]",
    ),
    "nhanes_seroprevalence": PrevalenceVariable(
        variable_type=VariableType.MEASUREMENT.value,
        percentage=0.54,
        number_of_participants=8153,
        country="United States",
        source="https://www.cdc.gov/nchs/nhanes/",
    ),
}