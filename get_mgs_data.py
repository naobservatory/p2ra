from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Callable, Generic, NewType, TypeVar

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


T = TypeVar("T")


@dataclass
class Tree(Generic[T]):
    data: T
    children: list[Tree[T]] = field(default_factory=list)

    def _helper(self, depth: int) -> str:
        _spacer = "."
        return f"{_spacer * depth}{self.data}\n" + "".join(
            c._helper(depth + 1) for c in self.children
        )

    def __str__(self) -> str:
        return self._helper(0)

    def __iter__(self):
        yield self
        for child in self.children:
            yield from child

    def __getitem__(self, val: T) -> Tree[T] | None:
        return get_subtree(self, val)


def get_subtree(tree: Tree[T], val: T) -> Tree[T] | None:
    """Depth-first search for taxid"""
    for subtree in tree:
        if subtree.data == val:
            return subtree
    else:
        return None


def tree_from_list(input: list) -> Tree:
    return Tree(data=input[0], children=[tree_from_list(c) for c in input[1:]])


# TODO: Test that these are inverse functions
def tree_to_list(tree: Tree) -> list:
    return list(tree.data, *(tree_to_list(c) for c in tree.children))


S = TypeVar("S")


# TODO: Test map rules
def map_tree(tree: Tree[T], f: Callable[[T], S]) -> Tree[S]:
    return Tree(f(tree.data), [map_tree(c, f) for c in tree.children])


def load_tax_tree(mgs_dir: Path) -> Tree[TaxID]:
    with open(mgs_dir / "dashboard/human_virus_tree.json") as data_file:
        data = json.load(data_file)
    return map_tree(tree_from_list(data), lambda x: TaxID(int(x)))


def make_count_tree(
    taxtree: Tree[TaxID], sample_counts: SampleCounts
) -> Tree[tuple[TaxID, Counter[Sample]]]:
    return map_tree(
        taxtree,
        lambda taxid: (taxid, Counter(sample_counts[taxid]))
        if taxid in sample_counts
        else (taxid, Counter()),
    )


def count_reads(
    taxtree: Tree[TaxID], sample_counts: SampleCounts
) -> Counter[Sample]:
    count_tree = make_count_tree(taxtree, sample_counts)
    return sum(
        (elem.data[1] for elem in count_tree),
        start=Counter(),
    )


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
