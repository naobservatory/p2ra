import json
import urllib.request
from collections import Counter
from dataclasses import dataclass, field
from datetime import date
from typing import NewType, Optional

from pydantic import BaseModel

from pathogen_properties import TaxID
from tree import Tree, tree_from_list

repo_url = "https://raw.githubusercontent.com/naobservatory/mgs-pipeline/main/"


@dataclass
class GitHubRepo:
    user: str
    repo: str
    branch: str
    url: str = field(init=False)

    def __post_init__(self):
        self.url = (
            f"https://raw.githubusercontent.com/"
            f"{self.user}/{self.repo}/{self.branch}/"
        )

    def get_file(self, path: str) -> str:
        file_url = self.url + path
        with urllib.request.urlopen(file_url) as response:
            if response.status == 200:
                return response.read()
            else:
                raise ValueError(
                    f"Failed to download {file_url}. "
                    f"Response status code: {response.status}"
                )


BioProject = NewType("BioProject", str)
Sample = NewType("Sample", str)


def load_bioprojects(repo: GitHubRepo) -> dict[BioProject, list[Sample]]:
    data = json.loads(repo.get_file("dashboard/metadata_bioprojects.json"))
    return {
        BioProject(bp): [Sample(s) for s in samples]
        for bp, samples in data.items()
    }


class SampleAttributes(BaseModel):
    country: str
    location: str
    fine_location: Optional[str] = None
    # Fixme: Not all the dates are real dates
    date: date | str
    reads: int


def load_sample_attributes(repo: GitHubRepo) -> dict[Sample, SampleAttributes]:
    data = json.loads(repo.get_file("dashboard/metadata_samples.json"))
    return {
        Sample(s): SampleAttributes(**attribs) for s, attribs in data.items()
    }


SampleCounts = dict[TaxID, dict[Sample, int]]


def load_sample_counts(repo: GitHubRepo) -> SampleCounts:
    data: dict[str, dict[str, int]] = json.loads(
        repo.get_file("dashboard/human_virus_sample_counts.json")
    )
    return {
        TaxID(int(taxid)): {Sample(sample): n for sample, n in counts.items()}
        for taxid, counts in data.items()
    }


def load_tax_tree(repo: GitHubRepo) -> Tree[TaxID]:
    data = json.loads(repo.get_file("dashboard/human_virus_tree.json"))
    return tree_from_list(data).map(lambda x: TaxID(int(x)))


def make_count_tree(
    taxtree: Tree[TaxID], sample_counts: SampleCounts
) -> Tree[tuple[TaxID, Counter[Sample]]]:
    return taxtree.map(
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
