import json
from pathlib import Path
from typing import Any, NewType

from pathogen_properties import TaxID
from pathogens import pathogens

Sample = NewType("Sample", str)


def load_samples(mgs_dir: Path, bioproject: str) -> list[Sample]:
    with open(mgs_dir / "metadata_bioprojects.json") as projects_file:
        data = json.load(projects_file)
    return [Sample(s) for s in data[bioproject]]


SampleAttributes = dict[str, str]


def load_sample_attributes(
    mgs_dir: Path, samples: list[Sample]
) -> dict[Sample, SampleAttributes]:
    with open(mgs_dir / "metadata_samples.json") as samples_file:
        data = json.load(samples_file)
    return {s: data[s] for s in samples}


def load_sample_counts(mgs_dir: Path) -> dict[TaxID, dict[Sample, int]]:
    with open(mgs_dir / "human_virus_sample_counts.json") as data_file:
        data: dict[str, dict[str, int]] = json.load(data_file)
    return {
        TaxID(int(taxid)): {Sample(sample): n for sample, n in counts.items()}
        for taxid, counts in data.items()
    }


pathogen = "sars_cov_2"
taxid = pathogens[pathogen].pathogen_chars.taxid
print(taxid)

# Rothman
bioproject = "PRJNA729801"
data_dir = Path("../mgs-pipeline/dashboard/")
samples = load_samples(data_dir, bioproject)
sample_attribs = load_sample_attributes(data_dir, samples)
for k, v in sample_attribs.items():
    print(k)
    print(v)
counts = load_sample_counts(data_dir)
print(counts[taxid])
