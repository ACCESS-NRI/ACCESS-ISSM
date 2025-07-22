#!/bin/bash
# run_mismip_first.sh  – ACCESS‑ISSM (vk83) edition
#
# Runs the MISMIP workflow just far enough to generate the
# Glen–vs‑ESTAR comparison plots, skipping the later “enhanced”
# or full‑Stokes experiments.
#
#      Mesh_generation → Parameterization → Transient_Steadystate
#      (optional) Transient_extrude
#      GlenSSA / GlenMOLHO / GlenHO
#      ESTARSSA / ESTARHO
#      analyse   ← comparison happens here
#
# Gadi job directives ---------------------------------------------------------
#PBS -P au88
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=48:00:00
#PBS -l jobfs=200GB
#PBS -l wd
#PBS -N mismip_first
#PBS -j oe
#------------------------------------------------------------------------------

set -eu  # fail hard on the first error

# ── 1.  Modules / software stack ─────────────────────────────────────────────
module purge
module use /g/data/vk83/modules            # << NEW: expose ACCESS‑ISSM modules
module load access-issm            #    (or pick a specific version)
# The access‑issm module sets:
#   * ISSM_DIR            → /g/data/vk83/...
#   * PATH / PYTHONPATH   → the matching Python 3.x executable & libs
#   * LD_LIBRARY_PATH     → PETSc, MUMPS, HDF5, NetCDF, etc.

# ── 2.  ISSM installation root ───────────────────────────────────────────────
# The module already exports $ISSM_DIR, so the next line is optional and safe
# to leave out unless you want to override it.
# export ISSM_DIR=$ISSM_DIR

# ── 3.  Define the “early” organiser steps we wish to execute  ---------------
early_steps='[1,2,3,8,9,10,11,18,20,31]'

# ── 4.  Patch the Python driver on‑the‑fly & run it  -------------------------
driver=./mismip_driver.py
tmp=$(mktemp mismipXXXX.py)

cp "$driver" "$tmp"
sed -i -E "s/^steps\s*=.*/steps = $early_steps/" "$tmp"

echo "Running first‑stage MISMIP workflow …"
python "$tmp"

rm "$tmp"
echo "Done – comparison plot and NetCDF outputs are in the Models_* folders."

