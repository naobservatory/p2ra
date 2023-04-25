import os.path
from dataclasses import dataclass
from enum import Enum
from typing import Optional

# Enums, short for enumerations, are a data type in Python used to represent a set of named values,
# which are typically used to define a set of related constants with unique names.

# TODO: Potentially create enum for country, state, etc. to ensure consistency in naming


class NAType(Enum):
    DNA = "DNA"
    RNA = "RNA"


class Enveloped(Enum):
    ENVELOPED = "enveloped"
    NON_ENVELOPED = "non_enveloped"


class VariableType(Enum):
    MEASUREMENT = "measurement"
    EXTERNAL_ESTIMATE = "external_estimate"
    NAO_ESTIMATE = "nao_estimate"


@dataclass(kw_only=True)
class PathogenChars:
    na_type: NAType
    enveloped: Enveloped
    taxid: int


@dataclass(kw_only=True)
class Variable:
    """An external piece of data"""

    source: Optional[str] = None
    country: Optional[str] = None
    state: Optional[str] = None
    county: Optional[str] = None
    number_of_participants: Optional[int] = None
    confidence_interval: Optional[tuple[float, float]] = None
    methods: Optional[str] = None
    # Either supply date, or start_date and end_date.
    # Dates can be any of: YYYY, YYYY-MM, or YYYY-MM-DD.
    date: Optional[str] = None
    start_date: Optional[str] = None
    end_date: Optional[str] = None
    is_target: Optional[bool] = False

    # Remember to recursively consider each input's inputs if defined.
    inputs: Optional[list["Variable"]] = None

    def _location(self):
        bits = []
        if self.county:
            bits.append(self.county)
        if self.state:
            bits.append(self.state)
        if self.country:
            bits.append(self.country)
        return ", ".join(bits)

    def _collect_locations(self, all_locations):
        location = self._location()
        if location:
            all_locations.add(location)
        if self.inputs and not self.is_target:
            for variable in self.inputs:
                variable._collect_locations(all_locations)

    def summarize_location(self, all_locations=None):
        all_locations = set()
        self._collect_locations(all_locations)
        return "; ".join(sorted(all_locations))

    def _collect_dates(self, all_dates):
        for date in [self.date, self.start_date, self.end_date]:
            if date:
                all_dates.add(date)
        if self.inputs and not self.is_target:
            for variable in self.inputs:
                variable._collect_dates(all_dates)

    def summarize_date(self):
        all_dates = set()
        self._collect_dates(all_dates)

        start_dates = set()
        end_dates = set()

        for date in all_dates:
            if len(date) == 4:
                start_dates.add("%s-01-01" % date)
                end_dates.add("%s-12-31" % date)
            elif len(date) == 7:
                start_dates.add("%s-01" % date)
                # Not technically correct, since some months are shorter, but
                # should be clear enough.
                end_dates.add("%s-31" % date)
            else:
                start_dates.add(date)
                end_dates.add(date)

        if not start_dates and not end_dates:
            return "no date"

        start_date = min(start_dates)
        end_date = max(end_dates)

        if start_date == end_date:
            return start_date

        start_year = start_date[:4]
        end_year = end_date[:4]

        if start_year != end_year:
            return "%s to %s" % (start_year, end_year)

        if start_date.endswith("-01-01") and end_date.endswith("-12-31"):
            return start_year

        start_month = start_date[5:7]
        end_month = start_date[5:7]

        if start_month != end_month:
            return start_year

        return "%s to %s" % (start_date, end_date)


@dataclass(kw_only=True)
class Population(Variable):
    """A number of people"""

    people: float
    # Make this specific enough that you won't accidentally pair it with the
    # wrong absolute prevalence or incidence.
    tag: str


@dataclass(kw_only=True)
class Scalar(Variable):
    scalar: float


@dataclass(kw_only=True)
class Prevalence(Variable):
    """What fraction of people have this pathogen at some moment"""

    infections_per_100k: float

    def scale(self, scalar: Scalar) -> "Prevalence":
        return Prevalence(
            infections_per_100k=self.infections_per_100k * scalar.scalar,
            inputs=[self, scalar],
        )

    def target(self, **kwargs) -> "Prevalence":
        return Prevalence(
            infections_per_100k=self.infections_per_100k,
            inputs=[self],
            is_target=True,
            **kwargs
        )


@dataclass(kw_only=True)
class PrevalenceAbsolute(Variable):
    """How many people had this pathogen at some moment"""

    infections: float
    # Make this specific enough that you won't accidentally pair it with the
    # wrong population.
    tag: str

    def to_rate(self, population: Population) -> Prevalence:
        assert self.tag == population.tag
        return Prevalence(
            infections_per_100k=self.infections * 100000 / population.people,
            inputs=[self, population],
        )


@dataclass(kw_only=True)
class SheddingDuration(Variable):
    days: float


@dataclass(kw_only=True)
class Number(Variable):
    """Generic number.  Use this for weird one-off things

    If some concept is being used by more than two estimates it should get a
    more specific Variable subclass."""

    number: float

    def per(self, other: "Number") -> Scalar:
        return Scalar(scalar=self.number / other.number, inputs=[self, other])


@dataclass(kw_only=True)
class IncidenceRate(Variable):
    """What fraction of people get this pathogen annually"""

    annual_infections_per_100k: float

    def to_prevalence(self, shedding_duration: SheddingDuration) -> Prevalence:
        return Prevalence(
            infections_per_100k=self.annual_infections_per_100k
            * shedding_duration.days
            / 365,
            inputs=[self, shedding_duration],
        )


@dataclass(kw_only=True)
class IncidenceAbsolute(Variable):
    """How many people get this pathogen annually"""

    annual_infections: float
    # Make this specific enough that you won't accidentally pair it with the
    # wrong population.
    tag: str

    def to_rate(self, population: Population) -> IncidenceRate:
        assert self.tag == population.tag
        return IncidenceRate(
            annual_infections_per_100k=self.annual_infections
            * 100000
            / population.people,
            inputs=[self, population],
        )


def prevalence_data_filename(filename):
    return os.path.join(os.path.dirname(__file__), "prevalence-data", filename)
