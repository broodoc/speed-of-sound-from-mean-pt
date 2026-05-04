#!/bin/bash

# Collect file sizes and paths
files=()
while IFS= read -r filepath; do
    filepath=$(echo "$filepath" | xargs)  # trim
    [[ ! -f "$filepath" ]] && continue
    size=$(stat -f "%z" "$filepath") || continue
    files+=("$size $filepath")
done < all_files.txt

# Sort by size (descending)
IFS=$'\n' sorted_files=($(printf "%s\n" "${files[@]}" | sort -nr))

# Init group accumulators
declare -a group_sizes
declare -a group_files
for i in {0..9}; do
    group_sizes[$i]=0
    group_files[$i]=""
done

# Greedy bin-packing: assign each file to group with min total size
for entry in "${sorted_files[@]}"; do
    size=${entry%% *}
    file=${entry#* }

    min_group=0
    for i in {1..9}; do
        if [ "${group_sizes[$i]}" -lt "${group_sizes[$min_group]}" ]; then
            min_group=$i
        fi
    done

    group_sizes[$min_group]=$(( group_sizes[$min_group] + size ))
    group_files[$min_group]+="$file"$'\n'
done

# Write each group to a group_X.txt file
for i in {0..9}; do
    echo -n "${group_files[$i]}" > "group_$((i+1)).txt"
done
