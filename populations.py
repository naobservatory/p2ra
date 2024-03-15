import datetime
import dataclasses
from typing import Optional

from pathogen_properties import Population, prevalence_data_filename

location_populations: list[tuple[str, dict[int, int]]] = []


def load_location_populations():
    if location_populations:
        return

    # Downloaded 2023-05-11 from
    # https://www2.census.gov/programs-surveys/popest/tables/2020-2022/counties/totals/co-est2022-pop.xlsx
    with open(prevalence_data_filename("Census-co-est2022-pop.tsv")) as inf:
        for line in inf:
            bits = line.strip().split("\t")
            if len(bits) != 5:
                continue

            location = bits[0]
            counts = {
                2020: int(bits[2].replace(",", "")),
                2021: int(bits[3].replace(",", "")),
                2022: int(bits[4].replace(",", "")),
            }
            location_populations.append((location, counts))


def us_population(
    year: int, county: Optional[str] = None, state: Optional[str] = None
) -> Population:

    if year in [2023, 2024]:
        pop_2022 = us_population(2022, county, state)
        return dataclasses.replace(
            pop_2022,
            parsed_start=datetime.date(year, 7, 1),
            parsed_end=datetime.date(year, 7, 1))
    
    if year not in [2020, 2021, 2022]:
        raise Exception("Unsupported year: %s" % year)
    
    total_people = 0
    # All estimates are July 1st, specifically.
    pop_date = "%s-07-01" % year
    source = "https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html"
    load_location_populations()
    for location, counts in location_populations:
        people = counts[year]
        if not county and not state and location == "United States":
            return Population(
                people=people,
                source=source,
                date=pop_date,
                country="United States",
            )

        if location == ".%s, %s" % (county, state):
            return Population(
                people=people,
                source=source,
                date=pop_date,
                country="United States",
                state=state,
                county=county,
            )

        if not county and location.endswith(", %s" % state):
            total_people += people

    if total_people == 0:
        raise Exception("county=%r, state=%r not found" % (county, state))

    # The only case where we total up is to get state populations.
    return Population(
        people=total_people,
        source=source,
        date=pop_date,
        country="United States",
        state=state,
    )
