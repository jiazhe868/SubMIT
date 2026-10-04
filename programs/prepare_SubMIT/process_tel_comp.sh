#!/bin/bash

# Navigate to the data directory
dir=$1
data_dir=${dir}/data
echo $data_dir $dir

# Check if the data directory exists
if [ ! -d "$data_dir" ]; then
    echo "Data directory does not exist: $data_dir"
    exit 1
fi

# Idempotency guard: integration is done IN PLACE, so a second run would
# double-integrate the data. Refuse to run twice on the same directory.
marker="$data_dir/.tel_comp_done"
if [ -f "$marker" ]; then
    echo "SKIP: $data_dir already integrated (remove $marker to force)"
    exit 0
fi
touch "$marker"

# Process *.t files: Integrate to displacement and keep the filenames
for file in "$data_dir"/*.t; do
    if [ -f "$file" ]; then
        sac << EOF
r $file
int
hp c 0.005 n 4 p 2
w $file
q
EOF
    fi
done

# Process *.z files: Integrate to displacement and keep the filenames, then save a velocity version with "vel_" prefix
for file in "$data_dir"/*.z; do
    if [ -f "$file" ]; then
        # Save the velocity version with "vel_" prefix
        velocity_file="$data_dir/vel_$(basename "$file")"
        cp "$file" "$velocity_file"

        # Integrate the original *.z file to displacement
        sac << EOF
r $file
int
hp c 0.005 n 4 p 2
w $file
q
EOF
    fi
done

