#!/bin/bash

echo "Running with shell: $SHELL"
echo "Bash version: $BASH_VERSION"

echo -=-=-=-=-=-=-=-=-=-=-

echo Shuffling .root files and listing in all_files.txt 🃏🔀
if [[ ! -f "$LIST" ]]; then
  echo "Shuffling .root files and listing in $LIST 🃏🔀"
  ls cleaned*root | sort -R > "$LIST"
  echo "$LIST created"
fi

echo -=-=-=-=-=-=-=-=-=-=-

#!/usr/bin/env bash
set -euo pipefail

#B=${1:-200}                 # number of bootstrap replicas
B=${1:-20}                 # number of bootstrap replicas
OUTDIR=${2:-bootstrap}      # output dir
LIST=${3:-all_files.txt}    # list of cleaned ROOT files

mkdir -p "$OUTDIR"

# Read file list into array
#mapfile -t files < <(grep -v '^\s*$' "$LIST")
files=()
while IFS= read -r line; do
  [[ -z "$line" ]] && continue
  files+=("$line")
done < "$LIST"

N=${#files[@]}

if [[ $N -lt 2 ]]; then
  echo "Need at least 2 files in $LIST"
  exit 1
fi

echo "Using N=$N files; generating B=$B bootstrap replicas in $OUTDIR"

# For each bootstrap replica: sample N indices with replacement, hadd them
for b in $(seq 1 "$B"); do
  rep_list="$OUTDIR/rep_${b}.txt"
  rep_root="$OUTDIR/boot_${b}.root"

  : > "$rep_list"

  # sample N draws with replacement
  for i in $(seq 1 "$N"); do
    r=$(( RANDOM % N ))
    echo "${files[$r]}" >> "$rep_list"
  done

  # hadd with duplicates allowed
  hadd -f "$rep_root" $(cat "$rep_list") > "$OUTDIR/hadd_${b}.log" 2>&1

  # basic sanity: ensure the output isn't tiny
  if [[ ! -s "$rep_root" ]]; then
    echo "Replica $b failed: $rep_root not created"
    exit 1
  fi

  echo "Made $rep_root"
done
