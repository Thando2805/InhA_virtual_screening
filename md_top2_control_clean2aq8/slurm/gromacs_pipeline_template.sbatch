#!/bin/bash
#SBATCH --job-name=inha_md_pilot
#SBATCH --partition=CHANGE_ME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=inha_md_%j.out
#SBATCH --error=inha_md_%j.err

set -euo pipefail

module purge
module load CHANGE_ME_GROMACS_MODULE

# Edit these names after you build the system.
SYSTEM_PREFIX="complex"
MDP_DIR="../mdp"

# Required files expected in the submit directory:
# - ${SYSTEM_PREFIX}.gro
# - topol.top
# - index.ndx

gmx grompp -f "${MDP_DIR}/minim.mdp" -c "${SYSTEM_PREFIX}.gro" -p topol.top -n index.ndx -o em.tpr -maxwarn 1
gmx mdrun -deffnm em

gmx grompp -f "${MDP_DIR}/nvt.mdp" -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt

gmx grompp -f "${MDP_DIR}/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt

gmx grompp -f "${MDP_DIR}/md_10ns_test.mdp" -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_10ns_test.tpr -maxwarn 1
gmx mdrun -deffnm md_10ns_test
