#!/usr/bin/env python3
"""
combine_topology.py
===================
Combines a GROMACS protein topology (from pdb2gmx) with a ligand topology
(from ACPYPE/GAFF2) into a single complex ready for solvation and MD.

Usage:
    python3 combine_topology.py

Run from the folder containing:
    vhl_processed.gro
    topol.top
    peptide_lig_H.acpype/peptide_lig_H_GMX.gro
    peptide_lig_H.acpype/peptide_lig_H_GMX.itp

Outputs:
    complex.gro   -- combined coordinates (protein + ligand)
    topol.top     -- updated topology including ligand parameters
"""

import sys
import os
import re

# -----------------------------------------------------------------------
# SETTINGS -- adjust if your filenames differ
# -----------------------------------------------------------------------
PROT_GRO    = "vhl_processed.gro"
LIG_GRO     = "peptide_lig_H.acpype/peptide_lig_H_GMX.gro"
LIG_ITP     = "peptide_lig_H.acpype/peptide_lig_H_GMX.itp"
LIG_POSRE   = "peptide_lig_H.acpype/posre_peptide_lig_H.itp"
TOPOL       = "topol.top"
OUT_GRO     = "complex.gro"
# -----------------------------------------------------------------------


def check_files():
    missing = [f for f in [PROT_GRO, LIG_GRO, LIG_ITP, TOPOL] if not os.path.exists(f)]
    if missing:
        sys.exit(f"Missing files:\n" + "\n".join(f"  {f}" for f in missing))
    print("All required files found.")


def get_molecule_name(itp_path):
    """Extract molecule name from ACPYPE-generated .itp file."""
    with open(itp_path) as f:
        lines = f.readlines()
    in_mol = False
    for line in lines:
        stripped = line.strip()
        if re.match(r'\[\s*moleculetype\s*\]', stripped):
            in_mol = True
            continue
        if in_mol and stripped and not stripped.startswith(';'):
            mol_name = stripped.split()[0]
            print(f"Ligand molecule name: '{mol_name}'")
            return mol_name
    sys.exit("Could not find [ moleculetype ] in ligand .itp file.")


def combine_gro(prot_gro, lig_gro, out_gro):
    """Merge protein and ligand .gro coordinate files."""
    with open(prot_gro) as f:
        prot_lines = f.readlines()
    with open(lig_gro) as f:
        lig_lines = f.readlines()

    prot_atoms = int(prot_lines[1].strip())
    lig_atoms  = int(lig_lines[1].strip())
    total      = prot_atoms + lig_atoms

    with open(out_gro, 'w') as f:
        f.write("VHL-peptide complex\n")
        f.write(f"{total}\n")
        # Protein atom lines (skip title + count + box vector)
        for line in prot_lines[2:2 + prot_atoms]:
            f.write(line)
        # Ligand atom lines
        for line in lig_lines[2:2 + lig_atoms]:
            f.write(line)
        # Box vectors from protein
        f.write(prot_lines[-1])

    print(f"complex.gro: {prot_atoms} protein + {lig_atoms} ligand = {total} atoms total")
    return prot_atoms, lig_atoms


def update_topology(topol, lig_itp, lig_posre, mol_name):
    """
    Edit topol.top to:
    1. Include the ligand .itp after the force field include
    2. Include ligand position restraints (ifdef POSRES_LIG)
    3. Add the ligand to the [ molecules ] section
    """
    with open(topol) as f:
        content = f.read()

    # ---- 1. Find the forcefield #include line and insert ligand itp after it
    ff_pattern = re.compile(r'(#include\s+"[^"]*forcefield\.itp")', re.MULTILINE)
    match = ff_pattern.search(content)
    if not match:
        sys.exit("Could not find forcefield #include in topol.top. "
                 "Check the file manually.")

    lig_itp_include = (
        f'\n; Ligand topology (GAFF2 via ACPYPE)\n'
        f'#include "{lig_itp}"\n'
    )

    # Only add if not already present
    if lig_itp not in content:
        content = content[:match.end()] + lig_itp_include + content[match.end():]
        print(f"Added #include \"{lig_itp}\" to topol.top")
    else:
        print(f"Ligand itp already included in topol.top -- skipping")

    # ---- 2. Add position restraint block for ligand (optional but recommended)
    if os.path.exists(lig_posre) and "POSRES_LIG" not in content:
        posre_block = (
            f'\n; Ligand position restraints\n'
            f'#ifdef POSRES_LIG\n'
            f'#include "{lig_posre}"\n'
            f'#endif\n'
        )
        # Insert before [ system ] or at end of itp include block
        sys_pattern = re.compile(r'(\[\s*system\s*\])', re.MULTILINE)
        sys_match = sys_pattern.search(content)
        if sys_match:
            content = content[:sys_match.start()] + posre_block + content[sys_match.start():]
        print(f"Added POSRES_LIG block for {lig_posre}")

    # ---- 3. Add ligand to [ molecules ] section
    mol_pattern = re.compile(r'(\[\s*molecules\s*\].*?)(\Z)', re.DOTALL | re.MULTILINE)
    mol_match = mol_pattern.search(content)
    if not mol_match:
        sys.exit("Could not find [ molecules ] section in topol.top.")

    mol_section = mol_match.group(0)

    if mol_name not in mol_section:
        # Append ligand entry at the end of the molecules section
        lig_entry = f"{mol_name}             1\n"
        content = content.rstrip() + "\n" + lig_entry
        print(f"Added '{mol_name}  1' to [ molecules ] section")
    else:
        print(f"'{mol_name}' already in [ molecules ] -- skipping")

    # Write updated topology
    with open(topol, 'w') as f:
        f.write(content)

    print(f"topol.top updated successfully.")


def main():
    print("=" * 60)
    print("Combining VHL protein + peptide ligand topologies")
    print("=" * 60)

    check_files()

    # 1. Get ligand molecule name from the .itp file
    mol_name = get_molecule_name(LIG_ITP)

    # 2. Combine .gro files
    prot_atoms, lig_atoms = combine_gro(PROT_GRO, LIG_GRO, OUT_GRO)

    # 3. Update topol.top
    update_topology(TOPOL, LIG_ITP, LIG_POSRE, mol_name)

    print("\n" + "=" * 60)
    print("SUCCESS. Next steps:")
    print("=" * 60)
    print(f"""
1. Create the simulation box:
   gmx_mpi editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt cubic

2. Add solvent:
   gmx_mpi solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

3. Prepare ions:
   gmx_mpi grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
   echo "SOL" | gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top \\
       -pname NA -nname CL -neutral

4. Then run energy minimization, NVT, NPT, and production as before.

Key difference from the 1YCR run:
   During equilibration, add -DPOSRES_LIG to the mdp 'define' line to
   keep the ligand restrained while the protein equilibrates.
""")


if __name__ == "__main__":
    main()
