#!/usr/bin/env bash
# run_pipeline.sh -- POSIX equivalent of run_pipeline.ps1
# Usage:  ./run_pipeline.sh  [python-executable]
set -euo pipefail
PY="${1:-python3}"

for f in model_scope_v1.txt decisions_final.csv rules_v3_unfitted.txt \
         build_network.py make_network_v2.py make_rules.py \
         attractor_analysis.py run_screen.py; do
  [ -f "$f" ] || { echo "missing required input: $f" >&2; exit 1; }
done
[ -d run01/raw ] || echo "WARNING: run01/raw not found; OmniPath will be re-downloaded" >&2
"$PY" --version

echo; echo "=== STEP 1: build network from OmniPath ==="
"$PY" build_network.py --anchors KRAS MDM2 \
  --sinks CASP3 CASP7 PARP1 CDKN1A --nodes model_scope_v1.txt \
  --min-resources 2 --min-references 1 --min-curation-effort 5 \
  --terminal-sinks --decisions decisions_final.csv --cache run01/raw -d final

echo; echo "=== STEP 2: literature corrections ==="
"$PY" make_network_v2.py --in final/edges.csv --out .

echo; echo "=== STEP 3: validate rules against network ==="
"$PY" make_rules.py --edges network_v2_edges.csv --nodes network_v2_nodes.csv \
  --check rules_v3_unfitted.txt --fixed KRAS MDM2 MDM4 \
  --export-maboss model_v3.bnd --export-boolnet model_v3.bnet

echo; echo "=== STEP 4: exact attractor analysis ==="
"$PY" attractor_analysis.py --rules rules_v3_unfitted.txt \
  --clamp control:KRAS=1,MDM2=1,MDM4=1 \
  --clamp MDM2_degraded:KRAS=1,MDM2=0,MDM4=1 \
  --clamp KRAS_degraded:KRAS=0,MDM2=1,MDM4=1 \
  --clamp both_degraded:KRAS=0,MDM2=0,MDM4=1 \
  --samples 120 --steps 6000 --out attractor_report.md

echo; echo "=== STEP 5: 32-variant screen ==="
"$PY" run_screen.py --rules rules_v3_unfitted.txt \
  --out-json screen_results.json --out-md screen_report.md

echo; echo "=== STEP 6: checksums ==="
sha256sum network_v2_edges.csv network_v2_nodes.csv network_v2_changelog.csv \
  rules_v3_unfitted.txt attractor_report.md screen_report.md \
  screen_results.json decisions_final.csv model_scope_v1.txt | tee checksums.txt
echo; echo "PIPELINE COMPLETE"
