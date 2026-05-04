#!/bin/bash

# Count how many output.*.root files there are
files=(output.*.root)
total=${#files[@]}

cleaned_files=(cleaned_output.*.root)
cleaned_total=${#cleaned_files[@]}

if [[ $cleaned_total -eq $total && $total -gt 0 ]]; then
    echo "✅ All $total DSTs are already cleaned. Skipping cleaning step."
    exit 0
fi

start_time=$(date +%s)
last_eta="--:--"

# Loop with indexed counter
i=1
for file in "${files[@]}"; do
    loop_start=$(date +%s)

    # cleaning macro
    root -l -q -b "clean_runid_histos.C(\"$file\")" > /dev/null 2>&1

    loop_end=$(date +%s)

    # Only calculate ETA after at least 10 files processed
    if (( i >= 10 )); then
        elapsed=$((loop_end - start_time))
        avg_time=$((elapsed / i))
        remaining=$((avg_time * (total - i)))

        mins=$((remaining / 60))
        secs=$((remaining % 60))
        time_str=$(printf "%02d:%02d" "$mins" "$secs")

        # Update ETA only if not 00:00
        if [[ "$time_str" != "00:00" ]]; then
            last_eta="$time_str"
        fi
    fi

    echo -ne "🧹 Cleaning file $i of $total: $file | ⏳ ETA: ${last_eta}   \r"

    ((i++))
done

# Final line so we don't overwrite the last one
echo -e "\n✅ Finished cleaning $total files."
