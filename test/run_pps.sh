#!/bin/bash
#===============================================================================
#  run_pps.sh  — ACCESS-ISSM (vk83) edition
#===============================================================================
#
#  Purpose
#  -------
#  Run the **pre-processing** (PPS) phase of the MISMIP workflow on NCI Gadi.
#  This prepares the experiment geometry and parameter files by executing only
#  the first two organizer stages:
#
#       1. Mesh_generation
#       2. Parameterization
#
#  These produce the triangulated domain mesh and initialize all input fields
#  (bed, surface, mask, rheology, forcings) for later transient or steady-state
#  model steps. The resulting `Models_*` directories can be reused for the
#  main workflow driver (e.g., `run_mismip_first.sh`).
#
#  Output
#  ------
#  • Mesh (.exp / .nc) files in ./Models_<resolution>/
#  • Fully parameterized ISSM model (.py / .nc)
#
#  Notes
#  -----
#  - The ACCESS-ISSM module provides all runtime dependencies (PETSc, MUMPS,
#    NetCDF, HDF5, MPI, Python).
#  - Edit PBS resource requests below to match mesh size and turnaround needs.
#  - Recommended queue: “normal” for small–medium pre-processing jobs.
#
#===============================================================================
#  PBS Job Directives
#===============================================================================
#PBS -P au88               # Project code (change to your NCI allocation)
#PBS -q normal             # Queue name
#PBS -l ncpus=8            # CPU cores requested
#PBS -l mem=32GB           # Memory allocation
#PBS -l walltime=04:00:00  # Wall-clock limit
#PBS -l jobfs=50GB         # Local node scratch space
#PBS -l wd                 # Run from current working directory
#PBS -N mismip_pps         # Job name (appears in qstat / Gadi dashboard)
#PBS -j oe                 # Combine stdout and stderr
#===============================================================================

set -eu  # Exit immediately on error or unset variable

#──────────────────────────────────────────────────────────────────────────────
# 1.  Load ACCESS-ISSM software stack
#──────────────────────────────────────────────────────────────────────────────
module purge
module use /g/data/vk83/modules
module load access-issm
#  This module defines:
#     • ISSM_DIR            → installation root
#     • PATH, PYTHONPATH    → correct Python environment
#     • LD_LIBRARY_PATH     → PETSc / MUMPS / HDF5 / NetCDF libraries

#──────────────────────────────────────────────────────────────────────────────
# 2.  (Optional) override ISSM_DIR or tie to spack install
#──────────────────────────────────────────────────────────────────────────────
# The module already defines ISSM_DIR; overriding is optional.
# Uncomment below if you want to tie to a Spack-managed ISSM installation.
export ISSM_DIR=$(spack location -i issm)

# (Optional) Attach to persistent shell session if you launched one interactively.
# source ~/.persistent-sessions/bashrc &>/dev/null || true

#──────────────────────────────────────────────────────────────────────────────
# 3.  Restrict workflow to PPS steps only
#──────────────────────────────────────────────────────────────────────────────
early_steps='[1,2]'   # 1: Mesh_generation, 2: Parameterization

driver=./access-issm/examples/mismip/mismip_driver.py
tmp=$(mktemp mismip_ppsXXXX.py)

# Patch the Python driver to execute only steps 1 & 2
cp "$driver" "$tmp"
sed -i -E "s/^steps\\s*=.*/steps = $early_steps/" "$tmp"

echo "=================================================================="
echo " ACCESS-ISSM Pre-processing (PPS) for MISMIP"
echo "=================================================================="
echo " Working directory : $(pwd)"
echo " Using ISSM_DIR     : ${ISSM_DIR:-<not set>}"
echo " Steps to execute   : $early_steps"
echo "------------------------------------------------------------------"

python "$tmp"
rm -f "$tmp"

echo "------------------------------------------------------------------"
echo " PPS finished – pre-processed inputs stored under:"
echo "   $PBS_O_WORKDIR/Models_*/"
echo "------------------------------------------------------------------"
