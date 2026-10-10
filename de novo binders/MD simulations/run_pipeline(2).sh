#!/bin/bash
# =============================================================================
# run_pipeline.sh
#
# GROMACS 100 ns MD pipeline -- CHARMM36 / cubic box / TIP3P / Cl- only
# Protocol based on: doi:10.1186/s43088-026-00765-8 (100 ns instead of 300 ns)
#
# USAGE (run from the folder containing your PDB):
#   bash run_pipeline.sh PMI-MDM2.pdb
#
# Preprocessing runs here directly (<1 min), then em -> nvt -> npt -> md
# are submitted as a SLURM dependency chain. Log off after submission.
# =============================================================================
set -euo pipefail

PDB_FILE="${1:?Usage: bash run_pipeline.sh <input.pdb>}"

# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------

WATER_MODEL="tip3p"
BOX_TYPE="cubic"      # paper specifies cubic box
BOX_DIST="1.0"        # minimum solute-to-wall distance (nm)
ION_PNAME="NA"
ION_NNAME="CL"        # paper: only Cl- added; -neutral handles the decision

GMX_MODULES="gcc/9.4.0-eewq4j6 openmpi/4.1.2-akxtzzl gromacs/2022-47qrtrj"
PARTITION="gpu"
EQ_TIME="02:00:00"    # NVT/NPT are 100 ps each -- finish in minutes
PROD_TIME="24:00:00"  # per production chunk
PROD_CHUNKS=2         # 2 x 24h is a safe ceiling for 100 ns;
                       # GROMACS exits cleanly once nsteps is reached

if [ ! -f "$PDB_FILE" ]; then
    echo "ERROR: '$PDB_FILE' not found in $(pwd)." >&2
    exit 1
fi

module purge
module load $GMX_MODULES

# ---------------------------------------------------------------------------
#
# Three key differences from AMBER/OPLS setups:
#   1. rcoulomb = rvdw = 1.2 nm  (not 1.0)
#   2. vdwtype = Cut-off  +  vdw-modifier = Force-switch  +  rvdw-switch = 1.0
#   3. NO DispCorr -- CHARMM36 was parameterized without dispersion correction
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
; minim.mdp - 50,000 steps steepest descent, CHARMM36
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
; nvt.mdp - 100 ps NVT at 300 K, CHARMM36
define                = -DPOSRES
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
tc-grps               = Protein Non-Protein
tau_t                 = 0.1 0.1
ref_t                 = 300 300
constraints           = h-bonds
continuation          = no
gen_vel               = yes
gen_temp              = 300
gen_seed              = -1
EOF

cat > npt.mdp << 'EOF'
; npt.mdp - 100 ps NPT at 300 K / 1 bar, CHARMM36
define                = -DPOSRES
integrator            = md
nsteps                = 50000
dt                    = 0.002
cutoff-scheme         = Verlet
coulombtype           = PME
rcoulomb              = 1.2
vdwtype               = Cut-off
vdw-modifier          = Force-switch
rvdw-switch           = 1.0
rvdw                  = 1.2
tcoupl                = V-rescale
tc-grps               = Protein Non-Protein
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
; md.mdp - 100 ns production, CHARMM36
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
tc-grps               = Protein Non-Protein
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

echo "[1/3] CHARMM36 .mdp files written (rcoulomb=rvdw=1.2 nm, Force-switch vdW, no DispCorr)"

# ---------------------------------------------------------------------------
# PREPROCESSING (runs directly here, no GPU needed, takes <1 min)
# ---------------------------------------------------------------------------

# -ignh strips any existing H atoms and lets GROMACS re-add them correctly
# (important for AlphaFold outputs or PDB structures with non-standard H names)
gmx_mpi pdb2gmx \
    -f "$PDB_FILE" \
    -o processed.gro \
    -p topol.top \
    -water "$WATER_MODEL" \

gmx_mpi editconf \
    -f processed.gro \
    -o boxed.gro \
    -c \
    -d "$BOX_DIST" \
    -bt "$BOX_TYPE"

gmx_mpi solvate \
    -cp boxed.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p topol.top

gmx_mpi grompp \
    -f ions.mdp \
    -c solvated.gro \
    -p topol.top \
    -o ions.tpr

# SOL: select the solvent group for genion to replace with ions
echo "SOL" | gmx_mpi genion \
    -s ions.tpr \
    -o solv_ions.gro \
    -p topol.top \
    -pname "$ION_PNAME" \
    -nname "$ION_NNAME" \
    -neutral

echo "[2/3] Preprocessing complete -> solv_ions.gro"

# ---------------------------------------------------------------------------
# WRITE SLURM STAGE SCRIPTS
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

gmx_mpi grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx_mpi mdrun -deffnm em -nb gpu -pme cpu -pin on
# Note: -update gpu is intentionally omitted here; the GPU leapfrog update
# only applies to the md integrator, not steepest-descent minimization.
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

gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
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

gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx_mpi mdrun -deffnm npt -nb gpu -pme gpu -update gpu -pin on
EOF

cat > md.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_md
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --time=$PROD_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

# grompp only on the first chunk; subsequent chunks reuse md.tpr
if [ ! -f md.tpr ]; then
    gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
fi

# -maxh: saves a checkpoint and exits cleanly 30 min before SLURM's hard kill
# so the next chunk can pick up exactly where this one left off
if [ -f md.cpt ]; then
    gmx_mpi mdrun -deffnm md -cpi md.cpt -nb gpu -pme gpu -update gpu -pin on -maxh 23.5
else
    gmx_mpi mdrun -deffnm md -nb gpu -pme gpu -update gpu -pin on -maxh 23.5
fi
EOF

chmod +x em.slurm nvt.slurm npt.slurm md.slurm
echo "[3/3] SLURM stage scripts written: em.slurm nvt.slurm npt.slurm md.slurm"
