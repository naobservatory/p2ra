import csv
import pandas as pd
from scipy.stats import gmean  # type: ignore
from math import log
from collections import defaultdict

PERCENTILES = ["5%", "25%", "50%", "75%", "95%"]


def reads_df() -> pd.DataFrame:
    df = pd.read_csv("input.tsv", sep="\t")
    return df


def rothman_fits_data() -> pd.DataFrame:
    data = defaultdict(list)
    for p in PERCENTILES:
        data[f"{p}"] = []

    with open("fits_summary.tsv") as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            if row["location"] == "Overall":
                continue
            if row["study"] != "rothman":
                continue
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
        if virus not in target_viruses:
            continue
        virus_df = df[df["virus"] == virus]
        htp_df = virus_df[virus_df["location"] == "HTP"]

        non_htp_df = virus_df[virus_df["location"] != "HTP"]

        gmean_variance["virus"].append(virus)
        for quantile in PERCENTILES:
            non_htp_quantile_gm = (gmean(non_htp_df[quantile].dropna()),)
            htp_quantile = gmean(htp_df[quantile].dropna())
            variance = float(htp_quantile - non_htp_quantile_gm)

            gmean_variance[f"variance_{quantile}"].append(round(variance, 2))

    return pd.DataFrame(gmean_variance)


def start():
    df_fits = rothman_fits_data()

    variance_df = compute_geo_mean_ratio(df_fits)

    variance_df.to_csv("rothman_variance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    start()
