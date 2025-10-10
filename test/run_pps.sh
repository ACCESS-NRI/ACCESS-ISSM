#!/bin/bash
# run_pps.sh – Pre‑processing wrapper for ACCESS‑ISSM MISMIP

#PBS -P au88               # ← change to your project code
#PBS -q normal
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l walltime=04:00:00
#PBS -l jobfs=50GB
#PBS -l wd
#PBS -N mismip_pps
#PBS -j oe

set -eu

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

# export ISSM_DIR=$(spack location -i issm)

# 2. (Optional) tie into a persistent‑session if you started one
# source ~/.persistent-sessions/bashrc &>/dev/null || true

# 3. Restrict workflow to early steps: Mesh_generation + Parameterization
early_steps='[1,2]'

driver=./access-issm/examples/mismip/mismip_driver.py

tmp=$(mktemp mismip_ppsXXXX.py)
cp "$driver" "$tmp"
sed -i -E "s/^steps\s*=.*/steps = $early_steps/" "$tmp"

python "$tmp"
rm "$tmp"

echo "PPS finished – inputs in $PBS_O_WORKDIR/Models_*/"
