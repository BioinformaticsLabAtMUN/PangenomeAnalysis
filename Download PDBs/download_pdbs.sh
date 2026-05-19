#!/usr/bin/env bash
set -euo pipefail

# USAGE: ./download_pdbs_api.sh <metadata.tsv> <base_out_dir> [batch_size]
# batch_size is accepted for compatibility but not used (API is per-ID)

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
download_pdbs_api.sh — Download AlphaFold DB PDBs using the EBI API.

USAGE:
  download_pdbs_api.sh <metadata.tsv> <base_out_dir> [batch_size]

INPUT:
  metadata.tsv : tab-separated file with UniProt accessions in column 1.
                Header row allowed.

OUTPUT:
  <base_out_dir>/alphafold_pdbs/AF-<ACC>-F1-model_v<latest>.pdb
  <base_out_dir>/alphafold_pdbs/logs/download.log
  <base_out_dir>/uniprot_to_pdb_mapping.tsv

REQUIREMENTS:
  - curl
  - python
EOF
  exit 0
fi

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
  echo "Usage: $0 <metadata.tsv> <base_out_dir> [batch_size]" >&2
  exit 1
fi

METADATA="$1"
BASE_OUT="$2"

PDB_DIR="${BASE_OUT}/alphafold_pdbs"
LOG_DIR="${PDB_DIR}/logs"
LOG="${LOG_DIR}/download.log"
mkdir -p "$PDB_DIR" "$LOG_DIR"
: > "$LOG"

command -v curl >/dev/null 2>&1 || { echo "ERROR: curl not found" >&2; exit 127; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python not found" >&2; exit 127; }

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

IDS_FILE="${WORKDIR}/ids.txt"

# Extract IDs (skip header if it looks like one), dedupe
cut -f1 "$METADATA" \
  | awk 'NR==1 && tolower($0) ~ /(uniprot|accession|protein|id)/ {next} {print}' \
  | sed '/^[[:space:]]*$/d' \
  | sort -u > "$IDS_FILE"

echo "→ Total unique UniProt IDs: $(wc -l < "$IDS_FILE")" | tee -a "$LOG"

MAP="${BASE_OUT}/uniprot_to_pdb_mapping.tsv"
echo -e "accession\tpdb_file\tstatus" > "$MAP"

while read -r ACC; do
  [[ -z "${ACC// }" ]] && continue

  RESP="${WORKDIR}/${ACC}.json"

  # Fetch API response (retry a bit; include a UA just in case)
  if ! curl -fsSL --retry 3 --retry-delay 2 \
      -H "User-Agent: alphafold-pdb-downloader/1.0" \
      "https://alphafold.ebi.ac.uk/api/prediction/${ACC}" \
      -o "$RESP" >>"$LOG" 2>&1; then
    echo -e "${ACC}\t(not found)\tAPI_failed" >> "$MAP"
    echo "FAIL API $ACC" >> "$LOG"
    continue
  fi

  # Parse pdbUrl + latestVersion safely (never crash on invalid JSON)
  read -r PDBURL LATEST <<<"$(python3 - "$RESP" <<'PY'
import json,sys
path=sys.argv[1]
try:
    with open(path,'r',encoding='utf-8') as f:
        data=json.load(f)
except Exception:
    print("", "")
    sys.exit(0)

if isinstance(data,list) and data:
    x=data[0]
    print(x.get("pdbUrl",""), str(x.get("latestVersion","")))
else:
    print("", "")
PY
)"

  if [[ -z "$PDBURL" ]]; then
    # Log a snippet of what we got (helps debug HTML/empty/etc.)
    SNIP="$(head -c 200 "$RESP" | tr '\n' ' ' | tr '\r' ' ')"
    echo -e "${ACC}\t(not found)\tAPI_invalid_json_or_no_pdbUrl" >> "$MAP"
    echo "FAIL parse/no_pdbUrl $ACC | resp_snip=${SNIP}" >> "$LOG"
    continue
  fi

  OUT="${PDB_DIR}/AF-${ACC}-F1-model_v${LATEST}.pdb"

  if [[ -s "$OUT" ]]; then
    echo -e "${ACC}\t$(basename "$OUT")\tSuccess" >> "$MAP"
    echo "SKIP exists $ACC -> $(basename "$OUT")" >> "$LOG"
    continue
  fi

  if curl -fL --retry 3 --retry-delay 2 -o "$OUT" "$PDBURL" >>"$LOG" 2>&1; then
    echo -e "${ACC}\t$(basename "$OUT")\tSuccess" >> "$MAP"
    echo "OK $ACC -> $(basename "$OUT")" >> "$LOG"
  else
    rm -f "$OUT"
    echo -e "${ACC}\t(not found)\tDownload_failed" >> "$MAP"
    echo "FAIL download $ACC url=$PDBURL" >> "$LOG"
  fi

done < "$IDS_FILE"

echo "Download complete."
echo "• PDBs: $PDB_DIR"
echo "• Mapping: $MAP"
echo "• Log: $LOG"

