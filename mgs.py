import json
import urllib.request
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date
from enum import Enum
from typing import NewType, Optional

from pydantic import BaseModel

from pathogen_properties import TaxID
from tree import Tree

MGS_REPO_DEFAULTS = {
    "user": "naobservatory",
    "repo": "mgs-pipeline",
    "ref": "data-2023-06-02",
}

BioProject = NewType("BioProject", str)
Sample = NewType("Sample", str)


rna_bioprojects = {
    "crits_christoph": BioProject("PRJNA661613"),
    "rothman": BioProject("PRJNA729801"),
    "spurbeck": BioProject("PRJNA924011"),
}


@dataclass
class GitHubRepo:
    user: str
    repo: str
    ref: str

    def get_file(self, path: str) -> str:
        file_url = (
            f"https://raw.githubusercontent.com/"
            f"{self.user}/{self.repo}/{self.ref}/{path}"
        )
        with urllib.request.urlopen(file_url) as response:
            if response.status == 200:
                return response.read()
            else:
                raise ValueError(
                    f"Failed to download {file_url}. "
                    f"Response status code: {response.status}"
                )


def load_bioprojects(repo: GitHubRepo) -> dict[BioProject, list[Sample]]:
    data = json.loads(repo.get_file("dashboard/metadata_bioprojects.json"))
    return {
        BioProject(bp): [Sample(s) for s in samples]
        for bp, samples in data.items()
    }


class Enrichment(Enum):
    VIRAL = "viral"
    PANEL = "panel"


class SampleAttributes(BaseModel):
    country: str
    state: Optional[str] = None
    county: Optional[str] = None
    location: str
    fine_location: Optional[str] = None
    # Fixme: Not all the dates are real dates
    date: date | str
    reads: int
    enrichment: Optional[Enrichment] = None
    method: Optional[str] = None


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
    return Tree.tree_from_list(data).map(lambda x: TaxID(int(x)))


def make_count_tree(
    taxtree: Tree[TaxID], sample_counts: SampleCounts
) -> Tree[tuple[TaxID, Counter[Sample]]]:
    return taxtree.map(
        lambda taxid: (taxid, Counter(sample_counts[taxid]))
        if taxid in sample_counts
        else (taxid, Counter()),
    )


def count_reads(
    taxtree: Tree[TaxID] | None, sample_counts: SampleCounts
) -> Counter[Sample]:
    if taxtree is None:
        return Counter()
    count_tree = make_count_tree(taxtree, sample_counts)
    return sum(
        (elem.data[1] for elem in count_tree),
        start=Counter(),
    )


@dataclass
class MGSData:
    bioprojects: dict[BioProject, list[Sample]]
    sample_attrs: dict[Sample, SampleAttributes]
    read_counts: SampleCounts
    tax_tree: Tree[TaxID]

    @staticmethod
    def from_repo(
        user=MGS_REPO_DEFAULTS["user"],
        repo=MGS_REPO_DEFAULTS["repo"],
        ref=MGS_REPO_DEFAULTS["ref"],
    ):
        repo = GitHubRepo(user, repo, ref)
        return MGSData(
            bioprojects=load_bioprojects(repo),
            sample_attrs=load_sample_attributes(repo),
            read_counts=load_sample_counts(repo),
            tax_tree=load_tax_tree(repo),
        )

    def sample_attributes(
        self, bioproject: BioProject, enrichment: Optional[Enrichment] = None
    ) -> dict[Sample, SampleAttributes]:
        samples = {
            s: self.sample_attrs[s] for s in self.bioprojects[bioproject]
        }
        if enrichment:
            return {
                s: attrs
                for s, attrs in samples.items()
                if attrs.enrichment == enrichment
            }
        else:
            return samples

    def total_reads(self, bioproject: BioProject) -> dict[Sample, int]:
        return {
            s: self.sample_attrs[s].reads for s in self.bioprojects[bioproject]
        }

    def viral_reads(
        self, bioproject: BioProject, taxids: Iterable[TaxID]
    ) -> dict[Sample, int]:
        viral_counts_by_taxid = {
            taxid: count_reads(self.tax_tree[taxid], self.read_counts)
            for taxid in taxids
        }
        return {
            s: sum(viral_counts_by_taxid[taxid][s] for taxid in taxids)
            for s in self.bioprojects[bioproject]
        }
