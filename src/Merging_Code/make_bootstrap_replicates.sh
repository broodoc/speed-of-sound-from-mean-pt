#!/usr/bin/env bash
set -euo pipefail

echo "Running with shell: ${SHELL:-unknown}"
echo "Bash version: ${BASH_VERSION:-unknown}"
echo "-=-=-=-=-=-=-=-=-=-=-"

# ----------------------------
# Args / defaults
# ----------------------------
B=${1:-10}                    # number of bootstrap replicas to ensure
OUTDIR=${2:-bootstrap}         # output dir
LIST=${3:-all_files.txt}       # list of cleaned ROOT files
MIN_KB=${4:-500}               # minimum size threshold in KB
MIN_BYTES=$((MIN_KB * 1024))

mkdir -p "$OUTDIR"

# ----------------------------
# Ensure file list exists
# ----------------------------
if [[ ! -f "$LIST" ]]; then
  echo "Shuffling .root files and listing in $LIST 🃏🔀"
  ls cleaned*root | sort -R > "$LIST"
  echo "$LIST created"
fi

# ----------------------------
# Read file list into array (bash 3.2 compatible; no mapfile)
# ----------------------------
files=()
while IFS= read -r line; do
  line="${line#"${line%%[![:space:]]*}"}"   # ltrim
  line="${line%"${line##*[![:space:]]}"}"   # rtrim
  [[ -z "$line" ]] && continue
  [[ ! -f "$line" ]] && continue
  files+=("$line")
done < "$LIST"

N=${#files[@]}
if [[ $N -lt 2 ]]; then
  echo "Need at least 2 valid files in $LIST"
  exit 1
fi

echo "Using N=$N cleaned ROOT files"
echo "Ensuring B=$B bootstrap replicas in $OUTDIR with size >= ${MIN_KB} KB"
echo "-=-=-=-=-=-=-=-=-=-=-"

# ----------------------------
# Helper: get file size in bytes (macOS + linux)
# ----------------------------
get_size_bytes() {
  local f="$1"
  # macOS (BSD stat)
  if stat -f %z "$f" >/dev/null 2>&1; then
    stat -f %z "$f"
  else
    # Linux (GNU stat)
    stat -c %s "$f"
  fi
}

# ----------------------------
# Helper: (re)make one replica b
# ----------------------------
make_replica() {
  local b="$1"
  local rep_list="$OUTDIR/rep_${b}.txt"
  local rep_root="$OUTDIR/boot_${b}.root"
  local rep_log="$OUTDIR/hadd_${b}.log"

  : > "$rep_list"

  # sample N draws with replacement
  for i in $(seq 1 "$N"); do
    r=$(( RANDOM % N ))
    echo "${files[$r]}" >> "$rep_list"
  done

  # hadd with duplicates allowed
  hadd -f "$rep_root" $(cat "$rep_list") > "$rep_log" 2>&1 || true

  if [[ ! -f "$rep_root" ]]; then
    echo "❌ Replica $b failed: $rep_root not created (see $rep_log)"
    return 1
  fi

  local sz
  sz=$(get_size_bytes "$rep_root")
  if [[ "$sz" -lt "$MIN_BYTES" ]]; then
    echo "❌ Replica $b too small (${sz} bytes < ${MIN_BYTES}); see $rep_log"
    return 1
  fi

  echo "✅ Replica $b OK (${sz} bytes)"
  return 0
}

# ----------------------------
# Main: check existing replicas; only recreate missing/small ones
# ----------------------------
for b in $(seq 1 "$B"); do
  rep_root="$OUTDIR/boot_${b}.root"

  if [[ -f "$rep_root" ]]; then
    sz=$(get_size_bytes "$rep_root")
    if [[ "$sz" -ge "$MIN_BYTES" ]]; then
      echo "✅ Found $rep_root (${sz} bytes) — keeping"
      continue
    else
      echo "⚠️  Found $rep_root but too small (${sz} bytes < ${MIN_BYTES}) — recreating only this one"
    fi
  else
    echo "⚠️  Missing $rep_root — creating"
  fi

  # Try to recreate; if it still fails, retry a few times (rare hadd hiccups)
  ok=false
  for attempt in 1 2 3; do
    echo "→ Building replica $b (attempt $attempt/3)"
    if make_replica "$b"; then
      ok=true
      break
    fi
  done

  if [[ "$ok" != true ]]; then
    echo "❌ Gave up on replica $b after 3 attempts. Check:"
    echo "   $OUTDIR/hadd_${b}.log"
    exit 1
  fi
done

echo "-=-=-=-=-=-=-=-=-=-=-"
echo "Done: all $B bootstrap ROOT files exist in $OUTDIR and are >= ${MIN_KB} KB."
