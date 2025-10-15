#!/bin/bash
#===============================================================================
#  run_mismip_first.sh  — ACCESS-ISSM (vk83) edition
#===============================================================================
#
#  Purpose
#  -------
#  Launches the initial MISMIP (Marine Ice Sheet Model Intercomparison Project)
#  workflow on NCI Gadi using the ACCESS-ISSM module environment.
#  This script runs only the “core” MISMIP stages required to generate the
#  Glen–vs–ESTAR comparison plots (SSA / HO / MOLHO), omitting later
#  enhanced-rheology or full-Stokes experiments for brevity and turnaround.
#
#  Workflow summary
#  ----------------
#      1. Mesh_generation
#      2. Parameterization
#      3. Transient_Steadystate  → relaxation to steady geometry
#      4. (optional) Transient_extrude
#      5. GlenSSA / GlenMOLHO / GlenHO     → baseline flow regimes
#      6. ESTARSSA / ESTARHO               → anisotropic rheology variants
#      7. analyse                          → plots + NetCDF export
#
#  Output
#  ------
#  • Model checkpoints in  ./Models_<resolution>/
#  • Comparison plots (SSA vs MOLHO vs HO vs ESTAR)
#  • NetCDF archives for post-processing / plotting
#
#  Notes
#  -----
#  - The ACCESS-ISSM environment provides all dependencies (PETSc, MUMPS,
#    NetCDF4, HDF5, MPI, Python, etc.) via a single module load.
#  - Adjust PBS resources below (cpus, memory, walltime) to suit mesh size.
#  - Script is restart-safe: reruns overwrite only temporary driver copies.
#
#===============================================================================
#  PBS Job Directives  (NCI Gadi)
#===============================================================================
#PBS -P au88                   # Project code (charge account)
#PBS -q normal                 # Queue (use 'express' for quick tests)
#PBS -l ncpus=32               # Total CPUs (1 node × 32 cores typical)
#PBS -l mem=64GB               # Total memory
#PBS -l walltime=48:00:00      # Wall-clock time limit
#PBS -l jobfs=200GB            # Local scratch (temporary files)
#PBS -l wd                     # Run from current working directory
#PBS -N mismip_first           # Job name
#PBS -j oe                     # Combine stdout + stderr
#===============================================================================

set -eu  # Exit on first error or undefined variable

#──────────────────────────────────────────────────────────────────────────────
# 1.  Load ACCESS-ISSM software stack
#──────────────────────────────────────────────────────────────────────────────
module purge
module use /g/data/vk83/modules
module load access-issm
#  This module defines:
#    • ISSM_DIR            – installation root
#    • PATH, PYTHONPATH    – to ISSM’s Python executable/libraries
#    • LD_LIBRARY_PATH     – PETSc / MUMPS / HDF5 / NetCDF dependencies

#──────────────────────────────────────────────────────────────────────────────
# 2.  (Optional) Explicitly confirm or override ISSM_DIR
#──────────────────────────────────────────────────────────────────────────────
# export ISSM_DIR=/g/data/vk83/<username>/issm-install

#──────────────────────────────────────────────────────────────────────────────
# 3.  Define which “early” organiser steps to run
#──────────────────────────────────────────────────────────────────────────────
#  These correspond to the numbered blocks in mismip_driver.py
early_steps='[3,5,6,9,10,11,12,13,14,15,16,17,18,19]'

#──────────────────────────────────────────────────────────────────────────────
# 4.  Patch the Python driver on-the-fly and execute
#──────────────────────────────────────────────────────────────────────────────
driver=./mismip_driver.py
tmp=$(mktemp mismipXXXX.py)

cp "$driver" "$tmp"
sed -i -E "s/^steps\\s*=.*/steps = $early_steps/" "$tmp"

echo "=================================================================="
echo " Starting first-stage MISMIP workflow on Gadi ..."
echo " Working directory : $(pwd)"
echo " Using ISSM_DIR     : ${ISSM_DIR:-<not set>}"
echo " Steps to execute   : $early_steps"
echo "=================================================================="

python "$tmp"

rm -f "$tmp"
echo "------------------------------------------------------------------"
echo " Done – Glen/ESTAR comparison plots and NetCDF outputs saved under"
echo "   ./Models_* directories."
echo "------------------------------------------------------------------"
