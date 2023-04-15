from enum import Enum
# Enums, short for enumerations, are a data type in Python used to represent a set of named values,
# which are typically used to define a set of related constants with unique names.

#Defining an enumerator for a pathogen's nucleic acid type
class NAType(Enum):
	DNA = 'DNA'
	RNA = 'RNA'

#Defining an enumerator for a pathogen's enveloped or non-enveloped status
class EnvelopedNonEnveloped(Enum):
	enveloped = 'enveloped'
	non_enveloped = 'non_enveloped'

