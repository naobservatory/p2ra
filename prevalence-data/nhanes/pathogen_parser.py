# File to download and parse NHANES data for Hepatitis A, B, and C Virus, HIV, HSV-1, HSV-2, and HPV.

import subprocess

url = "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPA_J.XPT"

save_path = "hepa_j.xpt"

subprocess.run(["wget", "-O", save_path, url], check=True)

for year in ["2017-2018", "2015-2016"]:
    url_w_year = f"https://wwwn.cdc.gov/Nchs/Nhanes/{year}/"
    save_path = f"{year}"
    for leaf in ["HEPA_J", "HEPBD", "HEPB_J", "HEPC_J", "HEPE_J", "HIV_J"]:
        try:
            subprocess.run(
                ["wget", "-O", save_path, url_w_year + leaf], check=True
            )
        except subprocess.CalledProcessError as e:
            print(
                f"{leaf}.xpt is not available for the NHANES cycle {year}: {e}"
            )
