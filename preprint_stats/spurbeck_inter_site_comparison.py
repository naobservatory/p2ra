import csv
import pandas as pd
from scipy.stats import gmean
from math import log
from collections import defaultdict

PERCENTILES = ["5%", "25%", "50%", "75%", "95%"]


def reads_df() -> pd.DataFrame:
    df = pd.read_csv("input.tsv", sep="\t")
    return df


def spurbeck_fits_data() -> pd.DataFrame:
    data = {
        "predictor_type": [],
        "virus": [],
        "study": [],
        "location": [],
        "enriched": [],
    }
    for p in PERCENTILES:
        data[f"{p}"] = []

    with open("fits_summary.tsv") as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            if row["location"] == "Overall":
                continue

            if row["study"] != "spurbeck":
                continue

            if row["location"] in ["E", "F", "G", "H"]:
                data["enriched"].append(True)
            else:
                data["enriched"].append(False)

            data["predictor_type"].append(row["predictor_type"])
            data["virus"].append(row["tidy_name"])
            data["study"].append(row["study"])
            data["location"].append(row["location"])
            for p in PERCENTILES:
                data[f"{p}"].append(abs(log(float(row[f"{p}"]), 10)))

    df = pd.DataFrame.from_dict(data)
    return df


def compute_geo_mean_ratio(df: pd.DataFrame) -> pd.DataFrame:
    target_viruses = [
        "Norovirus (GI)",
        "Norovirus (GII)",
        "SARS-COV-2",
        "MCV",
        "JCV",
        "BKV",
    ]
    gmean_variance = defaultdict(list)
    for virus in df["virus"].unique():
        print(virus)
        if virus not in target_viruses:
            continue
        virus_df = df[df["virus"] == virus]
        enriched_virus_df = virus_df[virus_df["enriched"]]
        non_enriched_virus_df = virus_df[~virus_df["enriched"]]
        gmean_variance["virus"].append(virus)
        for quantile in PERCENTILES:
            enriched_gm = gmean(enriched_virus_df[quantile].dropna())
            non_enriched_gm = (
                gmean(non_enriched_virus_df[quantile].dropna()),
            )
            variance = float(enriched_gm - non_enriched_gm)
            print(variance)

            gmean_variance[f"variance_{quantile}"].append(round(variance, 2))

    return pd.DataFrame(gmean_variance)


def start():
    df_fits = spurbeck_fits_data()

    variance_df = compute_geo_mean_ratio(df_fits)

    variance_df.to_csv("spurbeck_variance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    start()
