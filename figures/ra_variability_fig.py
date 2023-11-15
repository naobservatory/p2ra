import csv
from dataclasses import dataclass
from math import log
import matplotlib.pyplot as plt  # type: ignore
import seaborn as sns

import numpy as np
import pandas as pd

PERCENTILES = [5, 25, 50, 75, 95]


def read_data() -> pd.DataFrame:
    data = {
        "predictor_type": [],
        "virus": [],
        "study": [],
        "location": [],
        "mean": [],
        "std": [],
        "min": [],
        "max": [],
    }
    for p in PERCENTILES:
        data[f"{p}%"] = []

    with open("fits_summary.tsv") as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            if row["location"] != "Overall":
                continue
            data["predictor_type"].append(row["predictor_type"])
            data["virus"].append(row["tidy_name"])
            data["study"].append(row["study"])
            data["location"].append(row["location"])
            data["mean"].append(float(row["mean"]))
            data["std"].append(float(row["std"]))
            data["min"].append(float(row["min"]))
            data["max"].append(float(row["max"]))
            for p in PERCENTILES:
                data[f"{p}%"].append(log(float(row[f"{p}%"]), 10))

    df = pd.DataFrame.from_dict(data)
    return df


def compute_diffs(df: pd.DataFrame) -> pd.DataFrame:
    viruses = df["virus"].unique()
    results_data = {
        "virus": [],
        "predictor_type": [],
        "diff_median": [],
        "low_diff": [],
        "high_diff": [],
    }

    for virus in viruses:
        virus_df = df[df["virus"] == virus]
        predictor_type = virus_df["predictor_type"].unique()

        min_median_index = virus_df["50%"].idxmin()
        max_median_index = virus_df["50%"].idxmax()

        diff_median = (
            virus_df.loc[min_median_index, "50%"]
            - virus_df.loc[max_median_index, "50%"]
        )
        low_diff = (
            virus_df.loc[min_median_index, "5%"]
            - virus_df.loc[max_median_index, "95%"]
        )
        high_diff = (
            virus_df.loc[min_median_index, "95%"]
            - virus_df.loc[max_median_index, "5%"]
        )
        results_data["virus"].append(virus)
        results_data["diff_median"].append(diff_median)
        results_data["low_diff"].append(low_diff)
        results_data["high_diff"].append(high_diff)
        results_data["predictor_type"].append(predictor_type)
    df = pd.DataFrame.from_dict(results_data)
    return df


df = compute_diffs(read_data())


def plot_custom_dot_plot(df):
    fig, ax = plt.subplots()
    # sort df by diff_median
    # df = df.sort_values(by=["diff_median"])
    scatter = ax.scatter(
        x=df["diff_median"],
        y=range(len(df)),
        alpha=0.6,
        edgecolors="w",
    )
    # add vertical lines across the plot.

    for i in range(len(df)):
        ax.plot(
            [df["low_diff"][i], df["high_diff"][i]],
            [i, i],
            color="k",
        )

    x_min, x_max = ax.get_xlim()
    for i in range(round(x_min), round(x_max)):
        ax.axvline(i, color="k", alpha=0.1)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["virus"])

    ax.set_title("RA(1%) Difference (median, 5%, 95%)")
    ax.set_xlabel("OOM difference between most variable study estimates")
    # Show the plot
    plt.show()
    plt.clf()


plot_custom_dot_plot(df)


# plot_dotwhiskers(df)
