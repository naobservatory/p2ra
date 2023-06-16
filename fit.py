#!/usr/bin/env python3
from pathlib import Path
from textwrap import dedent

import numpy as np
import pandas as pd

import stats
from mgs import MGSData, TaxID, rna_bioprojects
from pathogens import predictors_by_taxid


def geom_mean(x: pd.DataFrame) -> float:
    return np.exp(np.mean(np.log(x)))


def print_summary(
    pathogen: str,
    taxids: frozenset[TaxID],
    predictor: str,
    study: str,
    ra_at_1in1000: pd.DataFrame,
) -> None:
    title = f"{pathogen} {list(taxids)} relative abundance at 1 in 1000 {predictor} in {study.title()}"
    percentiles = [5, 25, 50, 75, 95]
    percentile_values = np.percentile(ra_at_1in1000, percentiles)
    d = 1
    sep = " " * 4
    output = f"""
    {"-" * len(title)}
    {title}
    {"-" * len(title)}
    Posterior arithmetic mean:
    {np.mean(ra_at_1in1000):.{d}e}
    Posterior geometric mean:
    {geom_mean(ra_at_1in1000):.{d}e}
    Posterior quantiles:
    {sep.join(f"{p:>{d+5}}%" for p in percentiles)}
    {sep.join(f"{x:.{d}e}" for x in percentile_values)}
    """
    print(dedent(output))


def start() -> None:
    outdir = Path("fits")
    outdir.mkdir(exist_ok=True)
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    mgs_data = MGSData.from_repo()
    for (
        pathogen_name,
        predictor_type,
        taxids,
        predictors,
    ) in predictors_by_taxid():
        for study, bioproject in rna_bioprojects.items():
            model = stats.build_model(
                mgs_data,
                bioproject,
                predictors,
                taxids,
                random_seed=1,
            )
            model.fit_model()
            df = model.dataframe
            assert df is not None
            ra_at_1in1000 = pd.pivot_table(
                df, index="draws", values=["ra_at_1in1000"]
            )
            print_summary(
                pathogen_name, taxids, predictor_type, study, ra_at_1in1000
            )

            prefix = f"{pathogen_name}-{list(taxids)}-{predictor_type}-{study}"
            fig_hist = model.plot_posterior_histograms()
            fig_hist.savefig(figdir / f"{prefix}-posthist.pdf")
            fig_viol = model.plot_violin()
            fig_viol.savefig(figdir / f"{prefix}-violin.pdf")
            fig_sigma_tau = model.plot_sigma_tau_density()
            fig_sigma_tau.savefig(figdir / f"{prefix}-sigma-tau.pdf")
            xys = [
                ("date", "viral_reads"),
                ("date", "predictor"),
                ("predictor", "viral_reads"),
            ]
            for x, y in xys:
                g = model.plot_posterior_samples(
                    x,
                    y,
                    style="county",
                    hue="fine_location",
                    hue_order=model.fine_locations,
                )
                if y == "predictor":
                    g.set(yscale="log")
                g.savefig(figdir / f"{prefix}-{y}-vs-{x}.pdf")

            df.to_csv(
                outdir / f"{prefix}.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )


if __name__ == "__main__":
    start()
