import json
from pathlib import Path
from typing import NewType

from pathogens import pathogens
from pathogen_properties import TaxID

for k, v in pathogens.items():
    print(k)
    print(v.pathogen_chars.taxid)

pathogen = "sars_cov_2"
taxid = pathogens[pathogen].pathogen_chars.taxid
print(taxid)

# Rothman
bioproject = "PRJNA729801"
data_dir = Path("../mgs-pipeline/dashboard/")

Sample = NewType("Sample", str)

with open(data_dir / "metadata_bioprojects.json") as projects_file:
    meta_project = json.load(projects_file)
samples = meta_project[bioproject]


with open(data_dir / "metadata_samples.json") as samples_file:
    meta_samples = json.load(samples_file)

sample_attributes = {s: meta_samples[s] for s in samples}


def load_sample_counts(mgs_dir: Path) -> dict[TaxID, dict[Sample, int]]:
    with open(mgs_dir / "human_virus_sample_counts.json") as data_file:
        data = json.load(data_file)
    return {
        TaxID(taxid): {Sample(sample): n for sample, n in counts.items()}
        for taxid, counts in data.items()
    }


counts = load_sample_counts(data_dir)
x: Sample
for k in counts[TaxID(5)].keys():
    x = k

print(counts)
