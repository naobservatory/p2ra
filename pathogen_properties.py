import calendar
import datetime
import os.path
import re
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


@dataclass(kw_only=True)
class PathogenChars:
    na_type: NAType
    enveloped: Enveloped
    taxid: TaxID


@dataclass(kw_only=True)
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

    # Remember to recursively consider each input's inputs if defined.
    inputs: Optional[list["Variable"]] = None

    def __post_init__(
        self,
        date: Optional[str],
        start_date: Optional[str],
        end_date: Optional[str],
    ):
        if date and (start_date or end_date):
            raise Exception("If you have start/end don't set date.")
        if (start_date and not end_date) or (end_date and not start_date):
            raise Exception("Start and end must go together.")
        if date:
            start_date = end_date = date

        if start_date:
            self.parsed_start = self._parse_date(start_date, "start")
        else:
            self.parsed_start = None

        if end_date:
            self.parsed_end = self._parse_date(end_date, "end")
        else:
            self.parsed_end = None

        if (
            self.parsed_start
            and self.parsed_end
            and self.parsed_start > self.parsed_end
        ):
            raise Exception("Start date can't be after end date")

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
                _, d = calendar.monthrange(int(y), int(m))

        return datetime.date(y, m, d)

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

    def _collect_dates(self, all_dates: set[datetime.date]):
        for date in [self.parsed_start, self.parsed_end]:
            if date:
                all_dates.add(date)
        if self.inputs and not self.is_target:
            for variable in self.inputs:
                variable._collect_dates(all_dates)

    def summarize_date(self) -> Optional[tuple[datetime.date, datetime.date]]:
        all_dates: set[datetime.date] = set()
        self._collect_dates(all_dates)

        if not all_dates:
            return None

        return min(all_dates), max(all_dates)


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

    active: Active
    infections_per_100k: float

    def scale(self, scalar: Scalar) -> "Prevalence":
        return Prevalence(
            infections_per_100k=self.infections_per_100k * scalar.scalar,
            inputs=[self, scalar],
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


@dataclass(kw_only=True)
class PrevalenceAbsolute(Variable):
    """How many people had this pathogen at some moment"""

    infections: float
    active: Active

    # Make this specific enough that you won't accidentally pair it with the
    # wrong population.
    tag: str

    def to_rate(self, population: Population) -> Prevalence:
        assert self.tag == population.tag
        return Prevalence(
            infections_per_100k=self.infections * 100000 / population.people,
            inputs=[self, population],
            active=self.active,
        )


@dataclass(kw_only=True)
class SheddingDuration(Variable):
    days: float


@dataclass(kw_only=True)
class YearsOfInfections(Variable):
    """Number to convert an incidence to a latent prevalence, if rather than the shedding duration we are interested in the number of years for which these infections have been occuring"""

    years: float


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

    # Any estimate derived from an incidence using a shedding duration must be an active estimate,
    # since multiplying by SheddingDuration calculates the amount of time the virus is actively shedding for,
    # which is not incorporated into a latent estimate.
    def to_prevalence(self, shedding_duration: SheddingDuration) -> Prevalence:
        return Prevalence(
            infections_per_100k=self.annual_infections_per_100k
            * shedding_duration.days
            / 365,
            inputs=[self, shedding_duration],
            active=Active.ACTIVE,
        )

    def to_latent_prevalence(self, years: YearsOfInfections) -> Prevalence:
        return Prevalence(
            infections_per_100k=self.annual_infections_per_100k * years.years,
            inputs=[self, years],
            active=Active.LATENT,
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
