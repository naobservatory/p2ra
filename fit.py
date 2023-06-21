#!/usr/bin/env python3
from pathlib import Path

import stats
from mgs import MGSData, rna_bioprojects
from pathogens import predictors_by_taxid


def start() -> None:
    outdir = Path("fits")
    outdir.mkdir(exist_ok=True)
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    mgs_data = MGSData.from_repo()
    with open("fit_summary.tsv", "w") as summary_file:
        print(
            *(
                ["pathogen", "taxids", "predictor_type", "study"]
                + stats.Model.summary_header()
            ),
            sep="\t",
            file=summary_file,
        )
        for (
            pathogen_name,
            predictor_type,
            taxids,
            predictors,
        ) in predictors_by_taxid():
            for study, bioproject in rna_bioprojects.items():
                prefix = (
                    f"{pathogen_name}-{list(taxids)}-{predictor_type}-{study}"
                )
                model = stats.build_model(
                    mgs_data,
                    bioproject,
                    predictors,
                    taxids,
                    random_seed=1,
                )
                model.fit_model()
                for line in model.summary():
                    print(
                        *(
                            [
                                pathogen_name,
                                list(taxids),
                                predictor_type,
                                study,
                            ]
                            + line
                        ),
                        sep="\t",
                        file=summary_file,
                        flush=True,
                    )
                model.plot_figures(figdir, prefix)
                assert model.dataframe is not None
                model.dataframe.to_csv(
                    outdir / f"{prefix}.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                )


if __name__ == "__main__":
    start()
