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


@dataclass
class PathogenChars:
    na_type: NAType
    enveloped: Enveloped
    taxid: int


@dataclass
class PrevalenceVariable:
    variable_type: VariableType
    percentage: float
    source: str
    country: str
    state: str = None
    county: str = None
    number_of_participants: int = None
    confidence_interval: Tuple[float, float] = None
    methods: str = None


@dataclass
class PrevalenceEstimator:
    value: float
    value_type: str
    source: str
    country: str
    year: int
    state: str = None
    county: str = None
