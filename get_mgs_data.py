from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import NewType

from pydantic import BaseModel

from pathogen_properties import TaxID
from pathogens import pathogens

Sample = NewType("Sample", str)


def load_samples(mgs_dir: Path, bioproject: str) -> list[Sample]:
    with open(
        mgs_dir / "dashboard/metadata_bioprojects.json"
    ) as projects_file:
        data = json.load(projects_file)
    return [Sample(s) for s in data[bioproject]]


class SampleAttributes(BaseModel):
    country: str
    location: str
    fine_location: str | None = None
    date: date
    reads: int


def load_sample_attributes(
    mgs_dir: Path, samples: list[Sample]
) -> dict[Sample, SampleAttributes]:
    with open(mgs_dir / "dashboard/metadata_samples.json") as samples_file:
        data = json.load(samples_file)
    return {s: SampleAttributes(**data[s]) for s in samples}


SampleCounts = dict[TaxID, dict[Sample, int]]


def load_sample_counts(mgs_dir: Path) -> SampleCounts:
    with open(
        mgs_dir / "dashboard/human_virus_sample_counts.json"
    ) as data_file:
        data: dict[str, dict[str, int]] = json.load(data_file)
    return {
        TaxID(int(taxid)): {Sample(sample): n for sample, n in counts.items()}
        for taxid, counts in data.items()
    }


@dataclass
class TaxTree:
    taxid: TaxID
    children: list[TaxTree] = field(default_factory=list)

    def _helper(self, depth: int) -> str:
        _spacer = "."
        return f"{_spacer * depth}{self.taxid}\n" + "".join(
            c._helper(depth + 1) for c in self.children
        )

    def __str__(self) -> str:
        return self._helper(0)

    def __iter__(self):
        yield self
        for child in self.children:
            yield from child

    def __getitem__(self, taxid: TaxID) -> TaxTree | None:
        return get_subtree(self, taxid)


def get_subtree(taxtree: TaxTree, taxid: TaxID) -> TaxTree | None:
    """Depth-first search for taxid"""
    for subtree in taxtree:
        if subtree.taxid == taxid:
            return subtree
    else:
        return None


def count_reads(
    taxtree: TaxTree, sample_counts: SampleCounts
) -> Counter[Sample]:
    if taxtree.taxid in sample_counts:
        c = Counter(sample_counts[taxtree.taxid])
    else:
        c = Counter()
    return sum(
        (count_reads(child, sample_counts) for child in taxtree.children),
        start=c,
    )


def _parse_taxtree(input: list) -> TaxTree:
    taxid = TaxID(int(input[0]))
    children = input[1:]
    return TaxTree(taxid=taxid, children=[_parse_taxtree(c) for c in children])


def load_tax_tree(mgs_dir: Path) -> TaxTree:
    with open(mgs_dir / "dashboard/human_virus_tree.json") as data_file:
        data = json.load(data_file)
    return _parse_taxtree(data)


bioproject = "PRJNA729801"  # Rothman
data_dir = Path("../mgs-pipeline/")
samples = load_samples(data_dir, bioproject)
sample_attribs = load_sample_attributes(data_dir, samples)
fine_locs = set(sample_attribs[s].fine_location for s in samples)
counts = load_sample_counts(data_dir)
taxtree = load_tax_tree(data_dir)

for pathogen in ["sars_cov_2", "norovirus"]:
    print(pathogen)
    taxid = pathogens[pathogen].pathogen_chars.taxid
    subtree = taxtree[taxid]
    if subtree:
        virus_counts = count_reads(subtree, counts)
        print("All", sum(virus_counts[s] for s in samples), sep="\t")
        for fine_loc in fine_locs:
            print(
                fine_loc,
                sum(
                    virus_counts[s]
                    for s in samples
                    if sample_attribs[s].fine_location == fine_loc
                ),
                sep="\t",
            )
