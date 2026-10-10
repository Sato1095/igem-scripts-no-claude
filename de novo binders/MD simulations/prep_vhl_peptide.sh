#!/bin/bash
# =============================================================================
# prep_vhl_peptide.sh -- Ligand preparation for the VHL / LA-Hyp-Y-Hle-P rerun
#
# Produces everything needed before run_pipeline.sh:
#     complex.gro   protein + peptide coordinates
#     topol.top     AMBER99SB-ILDN protein + GAFF2 peptide, combined
#
# WHY THE PEPTIDE IS TREATED AS A LIGAND
#   LA-Hyp-Y-Hle-P contains two non-canonical residues -- 4-hydroxyproline
#   and homoleucine. Homoleucine is absent from AMBER99SB-ILDN entirely, so
#   pdb2gmx cannot build a topology for the peptide as a polymer. It is
#   instead parametrised as a single small molecule with GAFF2 via ACPYPE,
#   and the two topologies are merged.
#
#   (4-hydroxyproline alone would not require this -- GROMACS ships an HYP
#   entry for AMBER99SB-ILDN. It is homoleucine that forces the ligand route.)
#
# INPUT
#   rank_1.pdb    Boltz-2 prediction of the VHL/peptide complex, with the
#                 peptide written as HETATM records labelled LIG
#
# Usage:
#   module load gcc/9.4.0-eewq4j6 openmpi/4.1.2-akxtzzl gromacs/2022-47qrtrj
#   module load amber/24
#   export AMBERHOME=/mnt/apps/prebuilt/amber/24/amber24
#   conda activate mdenv
#   bash prep_vhl_peptide.sh rank_1.pdb
# =============================================================================
set -euo pipefail

INPUT=${1:-rank_1.pdb}
GMX="gmx_mpi"
FF="6"        # AMBER99SB-ILDN
WATER="1"     # TIP3P

echo "============================================================"
echo "VHL / LA-Hyp-Y-Hle-P  --  ligand preparation"
echo "  input: $INPUT"
echo "============================================================"

# ── Checks ───────────────────────────────────────────────────────────────────
[ -f "$INPUT" ] || { echo "ERROR: $INPUT not found"; exit 1; }
: "${AMBERHOME:?AMBERHOME not set -- module load amber/24 first}"
command -v acpype >/dev/null || {
    echo "ERROR: acpype not found."
    echo "  pip cache purge && pip install --no-cache-dir acpype"
    exit 1; }
command -v obabel >/dev/null || {
    echo "ERROR: obabel not found (normally bundled with acpype)"; exit 1; }

# ── 1. Split protein and peptide ─────────────────────────────────────────────
# Boltz-2 writes the peptide as HETATM (chain B, residue LIG) and VHL as ATOM.
echo ""
echo "[1/6] Splitting protein and peptide"

grep "^ATOM"   "$INPUT" > vhl_protein.pdb; echo "END" >> vhl_protein.pdb
grep "^HETATM" "$INPUT" > peptide_lig.pdb; echo "END" >> peptide_lig.pdb

NP=$(grep -c '^ATOM'   vhl_protein.pdb || true)
NL=$(grep -c '^HETATM' peptide_lig.pdb || true)
echo "  VHL     : ${NP} atoms"
echo "  peptide : ${NL} atoms"
[ "$NL" -gt 0 ] || {
    echo "ERROR: no HETATM records found."
    echo "  The peptide may be written as ATOM in this prediction."
    echo "  Check with:  grep HETATM $INPUT | head"
    exit 1; }

# ── 2. Add hydrogens ─────────────────────────────────────────────────────────
# Boltz-2 output is heavy-atom only; ACPYPE needs explicit hydrogens to
# assign charges.
echo ""
echo "[2/6] Adding hydrogens to the peptide"
obabel peptide_lig.pdb -h -O peptide_lig_H.pdb 2>/dev/null
NH=$(grep -c '^ATOM\|^HETATM' peptide_lig_H.pdb || true)
echo "  peptide with H : ${NH} atoms"

# ── 3. GAFF2 parameters ──────────────────────────────────────────────────────
# -c gas : Gasteiger charges. Faster than AM1-BCC and adequate here, since
#          the peptide is neutral and the comparison is structural rather
#          than a free-energy calculation. Use -c bcc if charges matter.
# -n 0   : net charge. LA-Hyp-Y-Hle-P is neutral -- verify for your construct.
echo ""
echo "[3/6] Generating GAFF2 parameters (ACPYPE)"
acpype -i peptide_lig_H.pdb -c gas -n 0 -a gaff2

ACP="peptide_lig_H.acpype"
for f in "${ACP}/peptide_lig_H_GMX.itp" "${ACP}/peptide_lig_H_GMX.gro"; do
    [ -f "$f" ] || { echo "ERROR: ACPYPE did not produce $f"; exit 1; }
done
echo "  ${ACP}/peptide_lig_H_GMX.itp"
echo "  ${ACP}/peptide_lig_H_GMX.gro"
echo "  ${ACP}/posre_peptide_lig_H.itp"

# ── 4. Protein topology ──────────────────────────────────────────────────────
echo ""
echo "[4/6] Building VHL topology (pdb2gmx)"
echo -e "${FF}\n${WATER}" | $GMX pdb2gmx \
    -f vhl_protein.pdb \
    -o vhl_processed.gro \
    -p topol.top \
    -i posre.itp \
    -ignh

echo "  vhl_processed.gro  topol.top  posre.itp"

# ── 5. Combine ───────────────────────────────────────────────────────────────
echo ""
echo "[5/6] Combining topologies"
[ -f combine_topology.py ] || {
    echo "ERROR: combine_topology.py not found in this directory"; exit 1; }
python3 combine_topology.py

[ -f complex.gro ] || { echo "ERROR: complex.gro not produced"; exit 1; }

# ── 6. Verify ────────────────────────────────────────────────────────────────
echo ""
echo "[6/6] Verification"
echo "  complex.gro atoms : $(sed -n '2p' complex.gro | tr -d ' ')"
echo ""
echo "  topol.top includes:"
grep -n '#include' topol.top | sed 's/^/    /'
echo ""
echo "  [ molecules ] section:"
grep -A 10 '^\[ molecules \]' topol.top | sed 's/^/    /'

echo ""
echo "============================================================"
echo "Ready for run_pipeline.sh"
echo ""
echo "Two things to check before running:"
echo "  - the ligand molecule name in [ molecules ] matches the .itp"
echo "  - grompp will need -maxwarn 1 for the mixed force field"
echo "    (AMBER99SB-ILDN protein + GAFF2 ligand)"
echo "============================================================"
