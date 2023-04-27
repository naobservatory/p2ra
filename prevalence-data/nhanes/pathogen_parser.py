# File to download and parse NHANES data for Hepatitis A, B, and C Virus, HIV, HSV-1, HSV-2, and HPV.

import subprocess
import os

# Create a new directory


for year in ["2015-2016", "2017-2018"]:
    if not os.path.exists(year):
        os.mkdir(year)

    if year == "2015-2016" return "I" else return="J"


    for leaf in [
        f"HEPA_{end}",
        f"HEPB_{end}",
        f"HEPC_{end}",
        f"HEPE_{end}",
        f"HIV_{end}",
        f"DEMO_{end}",
        f"HSV_{end}",
        f"HPVSWR_{end}",
        f"HPVP_{end}",
        f"OHRPV_{end}",
    ]:
        save_path = f"{year}/{leaf}"
        url_w_year = f"https://wwwn.cdc.gov/Nchs/Nhanes/{year}/{leaf}.XPT"
        try:
            subprocess.run(["wget", "-O", save_path, url_w_year], check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"{leaf}.xpt is not available for the NHANES cycle {year}: {e}"
            )
