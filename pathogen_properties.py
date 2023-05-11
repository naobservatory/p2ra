import calendar
import datetime
import os.path
import re
from collections.abc import Iterable
from dataclasses import InitVar, dataclass, field
from enum import Enum
from typing import NewType, Optional

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


class Active(Enum):
    ACTIVE = "Active"
    LATENT = "Latent"


TaxID = NewType("TaxID", int)


@dataclass(kw_only=True, eq=True, frozen=True)
class PathogenChars:
    na_type: NAType
    enveloped: Enveloped
    # Set exactly one of taxid or taxids; read taxids.
    #
    # Normally you only should set taxid.  Set taxids in cases like the flu
    # where surveillance generally conflates Flu A and Flu B but they don't
    # form a clade.
    taxid: InitVar[Optional[TaxID]] = None
    taxids: Optional[frozenset[TaxID]] = None
    # If we produce any estimates more specific than the overall taxid,
    # subtaxids will contain all the secondary taxonomic ids we can generate.
    subtaxids: frozenset[TaxID] = frozenset()

    def __post_init__(self, taxid: Optional[TaxID]):
        assert bool(taxid) ^ bool(self.taxids)  # Exactly one should be set.
        if taxid:
            # A python wart is that frozen dataclasses don't have an exception
            # for __post_init__, and it thinks assignments here are mutation
            # instead of initialization.  That's why we're assigning with
            # __setattr__. See
            # https://stackoverflow.com/questions/53756788/how-to-set-the-value-of-dataclass-field-in-post-init-when-frozen-true/54119384#54119384
            object.__setattr__(self, "taxids", frozenset([taxid]))


def days_in_month(year: int, month: int) -> int:
    _, last_day = calendar.monthrange(year, month)
    return last_day


@dataclass(kw_only=True, eq=True, frozen=True)
class Variable:
    """An external piece of data"""

    source: Optional[str] = None
    country: Optional[str] = None
    state: Optional[str] = None
    county: Optional[str] = None
    number_of_participants: Optional[int] = None
    confidence_interval: Optional[tuple[float, float]] = None
    coverage_probability: Optional[float] = None
    methods: Optional[str] = None
    # Either supply date, or start_date and end_date.
    # Dates can be any of: YYYY, YYYY-MM, or YYYY-MM-DD.
    date: InitVar[Optional[str]] = None
    start_date: InitVar[Optional[str]] = None
    end_date: InitVar[Optional[str]] = None
    parsed_start: Optional[datetime.date] = field(init=False)
    parsed_end: Optional[datetime.date] = field(init=False)
    is_target: Optional[bool] = False
    taxid: Optional[TaxID] = None
    inputs: InitVar[Optional[Iterable["Variable"]]] = None
    all_inputs: set["Variable"] = field(init=False)

    def __post_init__(
        self,
        date: Optional[str],
        start_date: Optional[str],
        end_date: Optional[str],
        inputs: Optional[Iterable["Variable"]],
    ):
        # See comment above about __post_init__ for why we're using __setattr__.
        if date and (start_date or end_date):
            raise Exception("If you have start/end don't set date.")
        if (start_date and not end_date) or (end_date and not start_date):
            raise Exception("Start and end must go together.")
        if date:
            start_date = end_date = date

        parsed_start = None
        if start_date:
            parsed_start = self._parse_date(start_date, "start")
        object.__setattr__(self, "parsed_start", parsed_start)

        parsed_end = None
        if end_date:
            parsed_end = self._parse_date(end_date, "end")
        object.__setattr__(self, "parsed_end", parsed_end)

        if (
            self.parsed_start
            and self.parsed_end
            and self.parsed_start > self.parsed_end
        ):
            raise Exception("Start date can't be after end date")

        all_inputs = set(inputs or [])
        for variable in set(inputs or []):
            all_inputs |= variable.all_inputs
        object.__setattr__(self, "all_inputs", frozenset(all_inputs))

    def _parse_date(self, date: str, start_or_end: str) -> datetime.date:
        y, m, d = None, None, None
        if y_match := re.findall("^(\d\d\d\d)$", date):
            (y,) = y_match
        elif ym_match := re.findall("^(\d\d\d\d)-(\d\d)$", date):
            ((y, m),) = ym_match
        elif ymd_match := re.findall("^(\d\d\d\d)-(\d\d)-(\d\d)$", date):
            ((y, m, d),) = ymd_match
        else:
            raise Exception("Unknown date format %s" % date)

        y = int(y)

        if m:
            m = int(m)
        else:
            m = {"start": 1, "end": 12}[start_or_end]

        if d:
            d = int(d)
        else:
            if start_or_end == "start":
                d = 1
            else:
                d = days_in_month(int(y), int(m))

        return datetime.date(y, m, d)

    # Returns country, state, county, or raises an error if there are
    # conflicting locations.  If you hit an error here you probably need a
    # target() call.
    def target_location(
        self,
    ) -> tuple[Optional[str], Optional[str], Optional[str]]:
        if self.is_target and self.country:
            return self.country, self.state, self.county

        inputs = self.all_inputs | set([self])

        countries = set(i.country for i in inputs if i.country)
        states = set(i.state for i in inputs if i.state)
        counties = set(i.county for i in inputs if i.county)

        country = state = county = None

        (country,) = countries
        if states:
            (state,) = states
        if counties:
            (county,) = counties

        return country, state, county

    def summarize_location(self) -> str:
        country, state, county = self.target_location()
        return ", ".join(x for x in [county, state, country] if x)

    def summarize_date(self) -> Optional[tuple[datetime.date, datetime.date]]:
        if self.is_target and self.parsed_start and self.parsed_end:
            return self.parsed_start, self.parsed_end

        try:
            return min(
                i.parsed_start
                for i in self.all_inputs | set([self])
                if i.parsed_start
            ), max(
                i.parsed_end
                for i in self.all_inputs | set([self])
                if i.parsed_end
            )
        except ValueError:
            return None


@dataclass(kw_only=True, eq=True, frozen=True)
class Taggable(Variable):
    # In cases where the location and date isn't enough to identify the
    # population, you can set a more specific tag to reduce errors.  For
    # example, tag="18-49yo".
    tag: Optional[str] = None

    def assert_comparable(self, other: "Taggable"):
        v1 = self
        v2 = other

        # While dates are optional in general, they're required for taggables.
        assert v1.parsed_start
        assert v2.parsed_start
        assert v1.parsed_end
        assert v2.parsed_end

        assert v1.country == v2.country
        assert v1.state == v2.state
        assert v1.county == v2.county

        # Normally everything has to match, but it's ok if one of them
        # has a more specific date as long as it's within a year; populations
        # don't change quickly.
        assert v1.parsed_start.year == v2.parsed_start.year
        assert v1.parsed_end.year == v2.parsed_end.year

        assert v1.tag == v2.tag


@dataclass(kw_only=True, eq=True, frozen=True)
class Population(Taggable):
    """A number of people"""

    people: float


@dataclass(kw_only=True, eq=True, frozen=True)
class Scalar(Variable):
    scalar: float


@dataclass(kw_only=True, eq=True, frozen=True)
class Prevalence(Variable):
    """What fraction of people have this pathogen at some moment"""

    active: Active
    infections_per_100k: float

    def __mul__(self, scalar: Scalar) -> "Prevalence":
        return Prevalence(
            infections_per_100k=self.infections_per_100k * scalar.scalar,
            inputs=[self, scalar],
            active=self.active,
        )

    def __add__(self: "Prevalence", other: "Prevalence") -> "Prevalence":
        assert self.active == other.active
        return Prevalence(
            infections_per_100k=self.infections_per_100k
            + other.infections_per_100k,
            inputs=[self, other],
            active=self.active,
        )

    def target(self, **kwargs) -> "Prevalence":
        return Prevalence(
            infections_per_100k=self.infections_per_100k,
            inputs=[self],
            is_target=True,
            active=self.active,
            **kwargs,
        )


@dataclass(kw_only=True, eq=True, frozen=True)
class PrevalenceAbsolute(Taggable):
    """How many people had this pathogen at some moment"""

    infections: float
    active: Active

    def to_rate(self, population: Population) -> Prevalence:
        self.assert_comparable(population)

        return Prevalence(
            infections_per_100k=self.infections * 100000 / population.people,
            inputs=[self, population],
            active=self.active,
        )


@dataclass(kw_only=True, eq=True, frozen=True)
class SheddingDuration(Variable):
    days: float


@dataclass(kw_only=True, eq=True, frozen=True)
class Number(Variable):
    """Generic number.  Use this for weird one-off things

    If some concept is being used by more than two estimates it should get a
    more specific Variable subclass."""

    number: float

    def __truediv__(self, other: "Number") -> Scalar:
        return Scalar(scalar=self.number / other.number, inputs=[self, other])


@dataclass(kw_only=True, eq=True, frozen=True)
class IncidenceRate(Variable):
    """What fraction of people get this pathogen annually"""

    annual_infections_per_100k: float

    # Any estimate derived from an incidence using a shedding duration must be
    # an active estimate, since multiplying by SheddingDuration calculates the
    # amount of time the virus is actively shedding for, which is not
    # incorporated into a latent estimate.
    def to_prevalence(self, shedding_duration: SheddingDuration) -> Prevalence:
        return Prevalence(
            infections_per_100k=self.annual_infections_per_100k
            * shedding_duration.days
            / 365,
            inputs=[self, shedding_duration],
            active=Active.ACTIVE,
        )


@dataclass(kw_only=True, eq=True, frozen=True)
class IncidenceAbsolute(Taggable):
    """How many people get this pathogen annually"""

    annual_infections: float

    def to_rate(self, population: Population) -> IncidenceRate:
        self.assert_comparable(population)

        return IncidenceRate(
            annual_infections_per_100k=self.annual_infections
            * 100000
            / population.people,
            inputs=[self, population],
        )


def prevalence_data_filename(filename):
    return os.path.join(os.path.dirname(__file__), "prevalence-data", filename)
