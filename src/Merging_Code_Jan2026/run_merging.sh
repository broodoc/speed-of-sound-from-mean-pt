#!/bin/bash

echo "Running with shell: $SHELL"
echo "Bash version: $BASH_VERSION"

echo -=-=-=-=-=-=-=-=-=-=-

echo Shuffling .root files and listing in all_files.txt 🃏🔀
ls cleaned_output*root | sort -R > all_files.txt
echo all_files.txt created

echo -=-=-=-=-=-=-=-=-=-=-

all_present_and_valid=true
total_size=0
existing_files=0
sizes=()

for i in $(seq 1 10); do
    file="merged_${i}.root"
    if [[ ! -f "$file" ]]; then
        echo "⚠️ Missing: $file"
        all_present_and_valid=false
        break
    fi

    size=$(stat -f %z "$file")
    if [[ "$size" -le 100 ]]; then
        echo "⚠️  $file is too small: ${size} bytes"
        all_present_and_valid=false
        break
    fi
    
    sizes[$i]=$size
    total_size=$((total_size + size))
    ((existing_files++))
done

if $all_present_and_valid; then
    avg_size=$((total_size / existing_files))
    min_allowed=$((avg_size / 2))

    for i in $(seq 1 10); do
        size=${sizes[$i]}
        if [[ "$size" -lt "$min_allowed" ]]; then
            echo "⚠️  $file is too small compared to average (${size} < ${min_allowed})"
            all_present_and_valid=false
            break
        fi
    done
fi

if $all_present_and_valid; then
    echo "✅ All merged_[1–10].root files exist and are valid. Skipping hadd step."
else
    echo "📦 Re-merging all groups..."

    echo Creating 10 groups of merged_x.txt files for partitioning 📂
    ./create_ten_groups.sh
    echo merged_x.txt files created

    echo -=-=-=-=-=-=-=-=-=-=-

    echo Merging files in each group. Creating merged_x.root files 🔗
    ./merge_ten.sh
    echo merged_x.root files created. Happy Bootstrapping! 👢
fi
