import os
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
    "repo": "mgs-restricted",
    "ref": "data-2023-07-21",
}

BioProject = NewType("BioProject", str)
Sample = NewType("Sample", str)


target_bioprojects = {
    "johnson": [BioProject("MJ-2024-02-12"),
                BioProject("MJ-2024-03-04")],
}


@dataclass
class GitHubRepo:
    user: str
    repo: str
    ref: str

    def get_file(self, path: str) -> str:
        with open(os.path.expanduser(
                f"~/code/{self.repo}/{path}")) as inf:
            return inf.read()


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
    location: Optional[str] = None
    fine_location: Optional[str] = None
    # Fixme: Not all the dates are real dates
    date: date | str
    reads: int
    edta_treated: Optional[bool|str] = False
    readlengths: Optional[dict] = None
    ribofrac: Optional[float] = None
    nuclease_treated: Optional[bool|str] = False
    sample_volume_ml: Optional[str] = None
    enrichment: Optional[Enrichment] = None
    concentration_method: Optional[str] = None
    method: Optional[str] = None
    airport: Optional[str] = None
    collection: Optional[str] = None



def load_sample_attributes(repo: GitHubRepo) -> dict[Sample, SampleAttributes]:
    data = json.loads(repo.get_file("dashboard/metadata_samples.json"))
    sa = {Sample(s): SampleAttributes(**attribs)
          for s, attribs in data.items()}
    for a in sa.values():
        a.location = "n/a"
        a.fine_location = "n/a"

    return sa


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
        lambda taxid: (
            (taxid, Counter(sample_counts[taxid]))
            if taxid in sample_counts
            else (taxid, Counter())
        ),
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
        return {
            s: attrs
            for s, attrs in samples.items()
            # remove them from the second run when it was True/False, but not
            # the first run when it was Yes/No.  The second run it was much to
            # agressive.
            if not attrs.edta_treated == True
        }

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
