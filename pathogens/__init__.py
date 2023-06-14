import importlib
import os

from pathogen_properties import NAType

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
