#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <folder_path>"
    exit 1
fi

folder_path="$1"
folder_name=$(basename "$folder_path")
raw_magnitude=$(echo "$folder_name" | awk -F'-' '{print substr($4, match($4, /[0-9]/))}')

case "$raw_magnitude" in
    ''|*[!0-9]*) echo "ERROR: cannot parse magnitude from folder name '$folder_name'" \
                      "(expected Wilber style ...-mwXY-...)" >&2; exit 1;;
esac

# Handle magnitude: add a decimal point if it is greater than 10
if [ "$raw_magnitude" -ge 10 ]; then
    magnitude=$(echo "$raw_magnitude" | awk '{print substr($0,1,1) "." substr($0,2)}')
else
    magnitude="$raw_magnitude"
fi

# List SAC files with .z extension
for sac_file in $(find "$folder_path/data/" -maxdepth 1 -name "*.z" | head -n 1); do
    if [ -f "$sac_file" ]; then
        # Extract header information using saclst
        sac_info=$(saclst kzdate kztime evlo evla evdp f "$sac_file")
        
        # Parse the information using awk
        kzdate=$(echo "$sac_info" | awk '{print $2}')
        kztime=$(echo "$sac_info" | awk '{print $3}')
        evlo=$(echo "$sac_info" | awk '{print $4}')
        evla=$(echo "$sac_info" | awk '{print $5}')
        evdp=$(echo "$sac_info" | awk '{print $6}')
        
        # Print the extracted information
        echo "$folder_name $kzdate $kztime $magnitude $evlo $evla $evdp"
    else
        echo "No .z SAC files found in $folder_path"
        exit 1
    fi
done

