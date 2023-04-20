from enum import Enum
from dataclasses import dataclass, asdict
from typing import Optional, Tuple

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
class PrevalenceVariable:
    variable_type: VariableType
    percentage: float
    source: str
    country: str
    state: Optional[str] = None
    county: Optional[str] = None
    number_of_participants: Optional[int] = None
    confidence_interval: Optional[Tuple[float, float]] = None
    methods: Optional[str] = None


@dataclass(kw_only=True)
class PrevalenceEstimator:
    value: float
    value_type: str
    source: str
    country: str
    year: int
    state: Optional[str] = None
    county: Optional[str] = None
