#!/bin/bash
# =============================================================================
# run_pipeline_lig.sh
#
# GROMACS 100 ns MD pipeline -- protein-ligand complex
#   AMBER99SB-ILDN protein + GAFF2 ligand / cubic box / TIP3P / Cl- only
#
# Same structure as run_pipeline.sh: preprocessing runs here, then
# em -> nvt -> npt -> md are submitted as a SLURM dependency chain.
#
# NONBONDED SETTINGS
#   rcoulomb = rvdw = 1.2 nm, vdw-modifier = Force-switch from 1.0 nm,
#   no analytical dispersion correction.
#
#   These are the same settings used for every other simulation in this
#   project and are retained here for consistency across the study. They
#   are not the AMBER reference values (1.0 nm plain cutoff with
#   DispCorr = EnerPres): the longer cutoff is more accurate in direct
#   space, force-switching improves energy conservation, and omitting
#   DispCorr slightly underestimates density. The effect on structural
#   observables -- RMSD, RMSF, Rg, contacts, buried surface area -- is
#   negligible, but absolute free energies would be affected, so any
#   MM-GBSA calculation on these trajectories should note the deviation.
#
#   Report these values as written, not the AMBER defaults.
#
# DIFFERENCE FROM run_pipeline.sh
#   Starts from complex.gro rather than a PDB. The peptide contains
#   homoleucine, which pdb2gmx cannot build, so the topology was assembled
#   separately (prep_vhl_peptide.sh) and must not be regenerated.
#
# USAGE (from the folder containing complex.gro and topol.top):
#   bash run_pipeline_lig.sh
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------
BOX_TYPE="cubic"
BOX_DIST="1.0"
ION_PNAME="NA"
ION_NNAME="CL"

GMX_MODULES="gcc/9.4.0-eewq4j6 openmpi/4.1.2-akxtzzl gromacs/2022-47qrtrj"
PARTITION="gpu"
EQ_TIME="02:00:00"
PROD_TIME="24:00:00"

# Mixing AMBER99SB-ILDN with GAFF2 raises a charge-sum warning on every grompp
MAXWARN=2

for f in complex.gro topol.top; do
    [ -f "$f" ] || { echo "ERROR: $f not found in $(pwd)" >&2; exit 1; }
done

module purge
module load $GMX_MODULES

# ---------------------------------------------------------------------------
# MDP FILES
# ---------------------------------------------------------------------------
cat > ions.mdp << 'EOF'
; ions.mdp
integrator        = steep
emtol             = 1000.0
emstep            = 0.01
nsteps            = 1000
nstlist           = 1
cutoff-scheme     = Verlet
coulombtype       = cutoff
rcoulomb          = 1.2
vdwtype           = Cut-off
vdw-modifier      = Force-switch
rvdw-switch       = 1.0
rvdw              = 1.2
pbc               = xyz
EOF

cat > minim.mdp << 'EOF'
; minim.mdp - 50,000 steps steepest descent
integrator        = steep
emtol             = 1000.0
emstep            = 0.01
nsteps            = 50000
nstlist           = 1
cutoff-scheme     = Verlet
coulombtype       = PME
rcoulomb          = 1.2
vdwtype           = Cut-off
vdw-modifier      = Force-switch
rvdw-switch       = 1.0
rvdw              = 1.2
pbc               = xyz
EOF

cat > nvt.mdp << 'EOF'
; nvt.mdp - 100 ps NVT at 300 K, protein and ligand both restrained
define                = -DPOSRES -DPOSRES_LIG
integrator            = md
nsteps                = 50000
dt                    = 0.002
nstxout-compressed    = 500
nstenergy             = 500
nstlog                = 500
cutoff-scheme         = Verlet
coulombtype           = PME
rcoulomb              = 1.2
vdwtype               = Cut-off
vdw-modifier          = Force-switch
rvdw-switch           = 1.0
rvdw                  = 1.2
tcoupl                = V-rescale
tc-grps               = Protein_LIG Water_and_ions
tau_t                 = 0.1 0.1
ref_t                 = 300 300
pcoupl                = no
constraints           = h-bonds
continuation          = no
gen_vel               = yes
gen_temp              = 300
gen_seed              = -1
EOF

cat > npt.mdp << 'EOF'
; npt.mdp - 100 ps NPT at 300 K / 1 bar
define                = -DPOSRES -DPOSRES_LIG
integrator            = md
nsteps                = 50000
dt                    = 0.002
nstxout-compressed    = 500
nstenergy             = 500
nstlog                = 500
cutoff-scheme         = Verlet
coulombtype           = PME
rcoulomb              = 1.2
vdwtype               = Cut-off
vdw-modifier          = Force-switch
rvdw-switch           = 1.0
rvdw                  = 1.2
tcoupl                = V-rescale
tc-grps               = Protein_LIG Water_and_ions
tau_t                 = 0.1 0.1
ref_t                 = 300 300
pcoupl                = Parrinello-Rahman
pcoupltype            = isotropic
tau_p                 = 2.0
ref_p                 = 1.0
compressibility       = 4.5e-5
constraints           = h-bonds
continuation          = yes
refcoord_scaling      = com
EOF

cat > md.mdp << 'EOF'
; md.mdp - 100 ns production, unrestrained
integrator            = md
dt                    = 0.002
nsteps                = 50000000   ; 100 ns
cutoff-scheme         = Verlet
coulombtype           = PME
rcoulomb              = 1.2
vdwtype               = Cut-off
vdw-modifier          = Force-switch
rvdw-switch           = 1.0
rvdw                  = 1.2
tcoupl                = V-rescale
tc-grps               = Protein_LIG Water_and_ions
tau_t                 = 0.1 0.1
ref_t                 = 300 300
pcoupl                = Parrinello-Rahman
pcoupltype            = isotropic
tau_p                 = 2.0
ref_p                 = 1.0
compressibility       = 4.5e-5
constraints           = h-bonds
continuation          = yes
nstxout-compressed    = 5000   ; 1 frame per 10 ps -> 10,000 frames total
nstenergy             = 5000
nstlog                = 5000
EOF

echo "[1/4] .mdp files written (rcoulomb=rvdw=1.2 nm, Force-switch vdW, no DispCorr)"

# ---------------------------------------------------------------------------
# PREPROCESSING -- no pdb2gmx; complex.gro and topol.top already exist
# ---------------------------------------------------------------------------
gmx_mpi editconf -f complex.gro -o boxed.gro -c -d "$BOX_DIST" -bt "$BOX_TYPE"

gmx_mpi solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

gmx_mpi grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr \
    -maxwarn $MAXWARN

echo "SOL" | gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top \
    -pname "$ION_PNAME" -nname "$ION_NNAME" -neutral

echo "[2/4] Preprocessing complete -> solv_ions.gro"

# ---------------------------------------------------------------------------
# INDEX -- GROMACS classifies a GAFF2 ligand as "Other", not "Protein", so
# the solute must be merged into one group for temperature coupling.
# Without this the ligand is thermostatted together with the solvent.
# ---------------------------------------------------------------------------
LIG=$(awk '/^\[ molecules \]/{f=1;next} f && $1!~/^;/ && NF>=2 &&
           $1!="SOL" && $1!="NA" && $1!="CL" && $1!~/^Protein/ {print $1; exit}' \
      topol.top)
[ -n "$LIG" ] || { echo "ERROR: could not identify ligand in [ molecules ]" >&2; exit 1; }

GRPS=$(printf 'q\n' | gmx_mpi make_ndx -f solv_ions.gro -o /tmp/probe_$$.ndx 2>&1 \
       | grep -E '^ *[0-9]+ ' || true)
rm -f /tmp/probe_$$.ndx
PROT_N=$(echo "$GRPS" | awk '$2=="Protein" {print $1; exit}')
LIG_N=$(echo  "$GRPS" | awk -v l="$LIG" '$2==l {print $1; exit}')
NEXT_N=$(echo "$GRPS" | wc -l)

[ -n "$LIG_N" ] || {
    echo "ERROR: ligand group '$LIG' not found. Groups were:" >&2
    echo "$GRPS" >&2; exit 1; }

printf '%s | %s\nname %s Protein_LIG\nq\n' "$PROT_N" "$LIG_N" "$NEXT_N" \
    | gmx_mpi make_ndx -f solv_ions.gro -o index.ndx > /dev/null 2>&1

echo "[3/4] index.ndx written -- ligand '$LIG' merged into Protein_LIG"

# ---------------------------------------------------------------------------
# SLURM STAGE SCRIPTS
# ---------------------------------------------------------------------------
cat > em.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_em
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --time=$EQ_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

gmx_mpi grompp -f minim.mdp -c solv_ions.gro -p topol.top -n index.ndx \\
    -o em.tpr -maxwarn $MAXWARN
gmx_mpi mdrun -deffnm em -nb gpu -pme cpu -pin on
# -update gpu omitted: the GPU leapfrog update applies to md, not steepest descent
EOF

cat > nvt.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_nvt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --time=$EQ_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx \\
    -o nvt.tpr -maxwarn $MAXWARN
gmx_mpi mdrun -deffnm nvt -nb gpu -pme gpu -update gpu -pin on
EOF

cat > npt.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_npt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --time=$EQ_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top \\
    -n index.ndx -o npt.tpr -maxwarn $MAXWARN
gmx_mpi mdrun -deffnm npt -nb gpu -pme gpu -update gpu -pin on
EOF

cat > md.slurm << EOF
#!/bin/bash
#SBATCH --partition=ampere
#SBATCH --job-name=gmx_md
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --time=10:00:00
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

# grompp only on the first chunk; subsequent chunks reuse md.tpr
if [ ! -f md.tpr ]; then
    gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top \\
        -n index.ndx -o md.tpr -maxwarn $MAXWARN
fi

# -maxh: checkpoint and exit cleanly 30 min before SLURM's hard kill
if [ -f md.cpt ]; then
    gmx_mpi mdrun -deffnm md -cpi md.cpt -nb gpu -pme gpu -update gpu \\
        -pin on -maxh 23.5
else
    gmx_mpi mdrun -deffnm md -nb gpu -pme gpu -update gpu -pin on -maxh 23.5
fi
EOF

chmod +x em.slurm nvt.slurm npt.slurm md.slurm
echo "[4/4] SLURM stage scripts written: em.slurm nvt.slurm npt.slurm md.slurm"
