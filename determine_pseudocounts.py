#!/usr/bin/env python3

import mgs
import pathogens

mgs_data = mgs.MGSData.from_repo()
import stats

for (
    pathogen_name,
    tidy_name,
    predictor_type,
    taxids,
    predictors,
) in pathogens.predictors_by_taxid():
    if not any(predictor.is_pseudocount for predictor in predictors):
        continue

    n_samples = 0
    n_pseudocount_samples = 0

    for study, bioprojects in mgs.target_bioprojects.items():
        for bioproject in bioprojects:
            enrichment = None if study == "brinch" else mgs.Enrichment.VIRAL
            chosen_predictors = {
                sample: stats.lookup_variables(attrs, predictors)
                for sample, attrs in mgs_data.sample_attributes(
                    bioproject, enrichment=enrichment
                ).items()
            }
            if all(ps == [] for ps in chosen_predictors.values()):
                continue
            for sample, preds in chosen_predictors.items():
                (predictor,) = preds
                n_samples += 1
                if predictor.is_pseudocount:
                    n_pseudocount_samples += 1

        print(
            tidy_name,
            study,
            n_pseudocount_samples,
            n_samples,
            "%.0f%%" % (100 * n_pseudocount_samples / n_samples),
        )
