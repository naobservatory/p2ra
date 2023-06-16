import importlib
import os
from typing import Generator

from mgs import TaxID
from pathogen_properties import Predictor, by_taxids

pathogens = {}
for pathogen_fname in os.listdir(os.path.dirname(__file__)):
    pathogen_name, ext = os.path.splitext(pathogen_fname)
    if pathogen_name == "__init__":
        continue
    if ext != ".py":
        continue

    pathogens[pathogen_name] = importlib.import_module(
        "pathogens.%s" % pathogen_name
    )


# Skip pathogens that don't have complete data yet.
# For Hepatitis B and C, we're waiting on extrapolating older estimates
# to the study period:
# https://github.com/naobservatory/p2ra/pull/154
skip = ["hbv", "hcv"]


def predictors_by_taxid() -> (
    Generator[tuple[str, str, frozenset[TaxID], list[Predictor]], None, None]
):
    pathogen_name: str
    predictor_type: str
    for pathogen_name, pathogen in pathogens.items():
        if pathogen_name in skip:
            continue
        for predictor_type, all_predictors in [
            ("incidence", pathogen.estimate_incidences()),
            ("prevalence", pathogen.estimate_prevalences()),
        ]:
            for taxids, predictors in by_taxids(
                pathogen.pathogen_chars,
                all_predictors,
            ).items():
                yield pathogen_name, predictor_type, taxids, predictors
