# HSV_1 prevalence calculation
# Written by Simon

# Importing enums and math
from enums import NAType, EnvelopedNonEnveloped
import math

# ============================
# Introduction and and Pathogen Characteristics
# ============================

background = """Herpes Simplex Virus 1 is a very common herpesvirus that causes oral herpes (CDC, https://www.
hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2), manifesting as cold sores or fever
blisters (https://www.hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2). After initial
infection and potential symptoms, most HSV-1 infections persist lifelong."""


pathogen_characteristics = {
	"na_type": NAType.dna,
	"enveloped_non_enveloped": EnvelopedNonEnveloped.non_enveloped,
	"taxid" : 10298
}

# ============================
# Variables for prevalence calculation
# ============================

prevalence_variables = {
	# In a the 2015-2016 NHANES survey, the seroprevalence of HSV-1 was 49.6% in sample of 3710 people:

	"seroprevalence": {
		"percentage": 0.496, #percentage prevalence
		"n": 3710,
		"confidence": 5, #subjective confidence from 1-6
		"source": "https://web.archive.org/web/20220707050306/https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSV_I.htm"
		},

	# In a 2015 study, 50 individuals self-collected tear and saliva samples twice daily which was tested for HSV 1
	# DNA. Individuals reported no prior ocular (eye) herpetic infection, but 23 reported prior herpetic disease.
	# The sample alse skewed toward African Americans (n=37; 78%). Across the 50-day study period, 49 out of 50
	# subjects shed HSV DNA atleast once.

	"tear_and_saliva_prevalence": {
		"percentage": 0.98, #percentage prevalence
		"n": 50,
		"confidence": 2, #subjective confidence from 1-6
		"source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1200985/#:~:text=Frequency%20of%20HSV%2D1%20DNA%20shedding%20over,swab%20data%20(subjects%2030%2C%2038%20not%20available)."
		}
}

# ============================
# Calculation of Aggregate Prevalence
# ============================

def calculate_aggregate_prevalence(prevalence_variables):
	# Store the variables for better readability
	sero_confidence = prevalence_variables["seroprevalence"]["confidence"]
	sero_percentage = prevalence_variables["seroprevalence"]["percentage"]
	tear_saliva_confidence = prevalence_variables["tear_and_saliva_prevalence"]["confidence"]
	tear_saliva_percentage = prevalence_variables["tear_and_saliva_prevalence"]["percentage"]

	# Calculate the geometric, weighted mean of the prevalences (https://en.wikipedia.org/wiki/Weighted_geometric_mean)

	# First, raise each prevalence to the power of its confidence
	seroprevalence_weighted = sero_percentage ** sero_confidence
	tear_saliva_prev_weighted = tear_saliva_percentage ** tear_saliva_confidence

	# Then multiply the two weighted prevalences
	prevalence_product = seroprevalence_weighted * tear_saliva_prev_weighted

	# Finally, take the nth root of the product, where n is the reciprocal of the number of weights.
	aggregate_prevalence = prevalence_product ** (1/7)

	return aggregate_prevalence


aggregate_prevalence = calculate_aggregate_prevalence(prevalence_variables)

# Making variables and functions available to other scripts
__all__ = ['pathogen_characteristics','prevalence_variables', 'aggregate_prevalence', 'background', 'calculate']