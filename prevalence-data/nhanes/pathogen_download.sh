#!/bin/bash

# Download xpt2csv.py
wget -O xpt2csv.py https://raw.githubusercontent.com/mirador/nhanes/master/xpt2csv.py

# Loop over years
for year in "2015-2016" "2017-2018"; do
    # Create year directory if it doesn't exist
    mkdir -p "$year"

    # Set the end variable based on the year
    if [ "$year" == "2017-2018" ]; then
        end="J"
    else
        end="I"
    fi

    # Loop over leaf names and download files
    for leaf in "HEPA_${end}" "HEPBD_${end}" "HEPB_S_${end}" "HEPC_${end}" "HEPE_${end}" "HIV_${end}" "DEMO_${end}" "HSV_${end}" "HPVSWR_${end}" "HPVP_${end}" "ORHPV_${end}"; do
        save_path="${year}/${leaf}.XPT"
        url_w_year="https://wwwn.cdc.gov/Nchs/Nhanes/${year}/${leaf}.XPT"

        # Download the file, skip if it already exists
        wget -nc -O "$save_path" "$url_w_year"
    done

    # Delete all files that contain "Oops" in their content
    grep -l "Oops" "${year}"/*.XPT | xargs rm

    # Run the xpt2csv.py script
    python xpt2csv.py "${year}" "${year}"

    # Remove all XPT files
    rm "${year}"/*.XPT
done

