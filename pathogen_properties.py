import calendar
import datetime
import math
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
    # In cases where an estimate is derived from multiple input variables with
    # different dates, set date_source to the Variable that represents the date
    # range this estimate is intended for.  For example, imagine we have:
    #
    #   population = Population(date="2020-04-01", ...)
    #   prevalence = AbsolutePrevalence(date="2020-08-01", ...)
    #   return prevalence.to_rate(population)
    #
    # There's no way for the consumer to know what date this estimate is for.
    # So instead we do:
    #
    #   return prevalence.to_rate(population, date_source=prevalence)
    #
    # parsed_start / parsed_end are set from date_source if supplied, otherwise
    # from start_date / end_date.
    date_source: InitVar[Optional["Variable"]] = None
    # Same deal as date_source: set this if otherwise it wouldn't be clear
    # what location an estimate would be for.
    location_source: InitVar[Optional["Variable"]] = None
    # Read these via get_dates(), which asserts that they're set.
    parsed_start: Optional[datetime.date] = None
    parsed_end: Optional[datetime.date] = None
    taxid: Optional[TaxID] = None
    inputs: InitVar[Optional[Iterable["Variable"]]] = None
    all_inputs: set["Variable"] = field(default_factory=set)

    def __post_init__(
        self,
        date: Optional[str],
        start_date: Optional[str],
        end_date: Optional[str],
        date_source: Optional["Variable"],
        location_source: Optional["Variable"],
        inputs: Optional[Iterable["Variable"]],
    ):
        # See comment above about __post_init__ for why we're using
        # __setattr__.
        if date and (start_date or end_date):
            raise Exception("If you have start/end don't set date.")
        if self.parsed_start and (date or start_date):
            raise Exception("Don't set both parsed_start and provide a date")
        if self.parsed_end and (date or end_date):
            raise Exception("Don't set both parsed_start and provide a date")
        if (start_date and not end_date) or (end_date and not start_date):
            raise Exception("Start and end must go together.")
        if date:
            start_date = end_date = date

        parsed_start = self.parsed_start
        if start_date:
            parsed_start = self._parse_date(start_date, "start")

        parsed_end = self.parsed_end
        if end_date:
            parsed_end = self._parse_date(end_date, "end")

        if date_source:
            assert date_source.parsed_start
            assert date_source.parsed_end
            parsed_start = date_source.parsed_start
            parsed_end = date_source.parsed_end

        if parsed_start and parsed_end and parsed_start > parsed_end:
            raise Exception("Start date can't be after end date")

        object.__setattr__(self, "parsed_start", parsed_start)
        object.__setattr__(self, "parsed_end", parsed_end)

        if location_source:
            object.__setattr__(self, "country", location_source.country)
            object.__setattr__(self, "state", location_source.state)
            object.__setattr__(self, "county", location_source.county)
        elif inputs:
            # If they didn't give us location information but there's location
            # information in our inputs, check that for consistency.  If it's
            # consistent use it, otherwise raise an error.

            countries = set(i.country for i in inputs if i.country)
            states = set(i.state for i in inputs if i.state)
            counties = set(i.county for i in inputs if i.county)

            country = self.country
            state = self.state
            county = self.county

            if not country:
                (country,) = countries
            if states and not state:
                (state,) = states
            if counties and not county:
                (county,) = counties

            object.__setattr__(self, "country", country)
            object.__setattr__(self, "state", state)
            object.__setattr__(self, "county", county)

        all_inputs = set(self.all_inputs or inputs or [])
        if date_source:
            all_inputs.add(date_source)
        if location_source:
            all_inputs.add(location_source)
        for variable in list(all_inputs):
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

    def get_dates(self) -> tuple[datetime.date, datetime.date]:
        assert self.parsed_start
        assert self.parsed_end
        return self.parsed_start, self.parsed_end

    def get_date(self) -> datetime.date:
        # Only call this on variables that you know represent a single-day
        # estimate.
        start, end = self.get_dates()
        assert start == end
        return start

    def get_location(
        self,
    ) -> tuple[Optional[str], Optional[str], Optional[str]]:
        return self.country, self.state, self.county

    def summarize_location(self) -> str:
        country, state, county = self.get_location()
        return ", ".join(x for x in [county, state, country] if x)


@dataclass(kw_only=True, eq=True, frozen=True)
class Taggable(Variable):
    # In cases where the location and date isn't enough to identify the
    # population, you can set a more specific tag to reduce errors.  For
    # example, tag="18-49yo".
    tag: Optional[str] = None

    def assert_comparable(self, other: "Taggable"):
        v1 = self
        v2 = other

        assert v1.country == v2.country
        assert v1.state == v2.state
        assert v1.county == v2.county

        # Normally everything has to match, but it's ok if one of them
        # has a more specific date as long as it's within a year; populations
        # don't change quickly.
        v1_start, v1_end = v1.get_dates()
        v2_start, v2_end = v2.get_dates()
        assert v1_start.year == v2_start.year
        assert v1_end.year == v2_end.year

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
            date_source=self,
            location_source=self,
        )

    def __add__(self: "Prevalence", other: "Prevalence") -> "Prevalence":
        assert self.active == other.active
        assert self.parsed_start == other.parsed_start
        assert self.parsed_end == other.parsed_end
        return Prevalence(
            infections_per_100k=self.infections_per_100k
            + other.infections_per_100k,
            inputs=[self, other],
            active=self.active,
            date_source=self,
            location_source=self,
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
            date_source=self,
            location_source=self,
        )


@dataclass(kw_only=True, eq=True, frozen=True)
class Number(Variable):
    """Generic number.  Use this for weird one-off things

    If some concept is being used by more than two estimates it should get a
    more specific Variable subclass."""

    number: float

    def __truediv__(self, other: "Number") -> Scalar:
        return Scalar(
            scalar=self.number / other.number,
            inputs=[self, other],
            date_source=self,
            location_source=self,
        )


@dataclass(kw_only=True, eq=True, frozen=True)
class IncidenceRate(Variable):
    """What fraction of people get this pathogen annually"""

    annual_infections_per_100k: float

    def __mul__(self, scalar: Scalar) -> "IncidenceRate":
        return IncidenceRate(
            annual_infections_per_100k=self.annual_infections_per_100k
            * scalar.scalar,
            inputs=[self, scalar],
            date_source=self,
            location_source=self,
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
            date_source=self,
            location_source=self,
        )


def prevalence_data_filename(filename):
    return os.path.join(os.path.dirname(__file__), "prevalence-data", filename)
