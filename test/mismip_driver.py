#!/usr/bin/env python3
"""
MISMIP end-to-end pipeline (Python/ISSM) — annotated
====================================================

This script mirrors a classic ISSM/MISMIP workflow (mesh → parameterize →
long relaxation → extrude → model-family experiments → analysis) and is laid
out in numbered **STEP** blocks you can toggle with an "organizer". The code is
heavily annotated to explain *what* happens and *why*, so you can treat it as a
living recipe for reproducing and extending your experiments on NCI Gadi or a
local machine.

Highlights
---------
- Switch between multiple model configurations (mesh resolution + friction law)
  by setting `modelnum` (1..8). The script derives a `modelname` string for
  bookkeeping and output paths.
- The **organizer** orchestrates persistence: each step can `savemodel`/`loadmodel`
  to/from a dedicated directory, enabling restartable workflows.
- Ready-to-use Gadi cluster block with sensible defaults (project, nodes, walltime),
  plus a single-node fallback for local runs.
- A set of optional modeling branches: SSA / HO / FS and "enhanced" rheologies
  (E or ESTAR variants) to compare flow approximations.

Usage
-----
1) Set environment and model choice
   - Ensure ISSM is installed and `ISSM_DIR` is set in the environment.
   - Choose which steps to run in `steps` and which variant via `modelnum`.
2) Run the script
   - `python mismip_pipeline_annotated.py`

Prerequisites
-------------
- ISSM Python interface on PATH (imports below must succeed)
- NetCDF4, NumPy, SciPy, Matplotlib available in the Python env
- Access to any datasets referenced in your `Mismip.py` parameterization

Notes
-----
- This script keeps the original structure and adds comments.
- Some lines (e.g., plotting/analysis) are placeholders and may need adaptation
  to your local ISSM build and data.
- For clarity, we *annotate* rather than refactor. Any obvious typos are noted
  inline as comments without changing behavior unless harmless.
"""

# ------------------------------
# Standard library & sanity checks
# ------------------------------
import os, sys
import math  # used for friction coefficients

# Ensure ISSM paths are correctly set (ISSM Python modules must be importable)
issm_dir = os.getenv('ISSM_DIR')
if not issm_dir:
    print("Error: ISSM_DIR environment variable is not set.")
    sys.exit(1)

# ------------------------------
# Third-party & ISSM imports
# ------------------------------
import numpy as np
from triangle import triangle
from model import *
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh
from bamg import bamg
from xy2ll import xy2ll
from plotmodel import plotmodel
from export_netCDF import export_netCDF
from loadmodel import loadmodel
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from socket import gethostname
from solve import solve
from ll2xy import ll2xy
from BamgTriangulate import BamgTriangulate
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from scipy.interpolate import griddata
from matplotlib import pyplot as plt
from gadi_spack import gadi
from organizer import organizer
from toolkits import toolkits
from bcgslbjacobioptions import bcgslbjacobioptions
from SetMOLHOBC import SetMOLHOBC

# ------------------------------
# 1) Choose which steps and model to run
# ------------------------------
# Define the list of step names to execute (must match org.perform(...) labels).
# Example single-step run:
# steps = ["Mesh_generation"]
# Example multi-step run:
# steps = ["Mesh_generation", "Parameterization", "Transient_Steadystate"]

# >>> EDIT ME <<<
# steps = [11]
# modelnum = 1

# ------------------------------
# 2) Derive a human-readable model name from the model number
# ------------------------------
# This controls directory/prefix names and helps keep outputs organized.
if 'modelnum' not in globals():
    raise ValueError("Please set 'modelnum' (1..8) near the top of the script.")

if modelnum == 1:
    modelname = '1km_viscous'
elif modelnum == 2:
    modelname = '2km_viscous'
elif modelnum == 3:
    modelname = '1km_coulomb'
elif modelnum == 4:
    modelname = '2km_coulomb'
elif modelnum == 5:
    modelname = '500m_viscous'
elif modelnum == 6:
    modelname = '500m_coulomb'
elif modelnum == 7:
    modelname = '200m_viscous'
elif modelnum == 8:
    modelname = '200m_coulomb'
else:
    raise ValueError('Invalid modelnum! Choose 1..8')

# ------------------------------
# Global run options
# ------------------------------
clustername = 'gadi'     # 'gadi' for NCI. Any other string triggers local fallback below.
printflag = True         # Reserved for any verbose printing you want to add

# ------------------------------
# Cluster configuration (Gadi vs local fallback)
# ------------------------------
# On Gadi we construct an ISSM-aware cluster object. Else we make a minimal
# single-node "generic" cluster (adjust np, mem, time as needed).
if clustername.lower() == 'gadi':

    cluster = gadi(
        'name',         'gadi.nci.org.au',
        'login',        'jh7060',                         # <<< EDIT: your NCI login
        'srcpath',      '/g/data/au88/jh7060/ISSM',       # <<< EDIT: where ISSM sources live
        'project',      'au88',                           # <<< EDIT: your NCI project code
        'numnodes',     1,
        'cpuspernode',  32,
        'time',         48*60,  # minutes
        'codepath',     '/g/data/au88/jh7060/spack/0.22/release/linux-rocky8-x86_64_v4/gcc-13.2.0/issm-4.24-v52nx3pfx7lpqwfldqlpspflk34wz756/bin',
        'executionpath','/scratch/au88/jh7060/issm_runs'
    )
    loadonly = 0  # run for real on the cluster
    lock     = 1  # use file locks for multi-job safety
else:
    # --- generic single-node/local fallback ---
    cluster        = generic()            # minimal cluster; no args
    cluster.name   = gethostname()        # label for logs
    cluster.np     = 32                   # number of MPI ranks
    cluster.mem    = 64                   # GB
    cluster.time   = 60                   # minutes
    cluster.queue  = 'local'
    lock     = 1
    loadonly = 1  # only prepare/run locally without submitting

# ------------------------------
# Organizer: centralizes persistence and step orchestration
# ------------------------------
# Models and results are saved under a model-specific directory with a prefix.
# The steps list controls which blocks below are executed.
if 'steps' not in globals():
    raise ValueError("Please set 'steps' near the top (list of step names to run).")

org = organizer(
     'repository',   './Models_'     + modelname,
     'prefix',       'mismip_'       + modelname + '_',
     'steps',        steps,
     'trunkprefix',  '34;47;2'  # free-form label; often used for provenance
)

# ============================================================================
#  STEP1: Mesh_generation
# ============================================================================
if org.perform('Mesh_generation'):
    # Start a blank ISSM model; then build an unstructured mesh from Domain.exp
    model = model()

    # Choose hmax based on modelnum (coarse → fine)
    if modelnum == 1 or modelnum == 3:
        md = bamg(model,
          'domain',       './Domain.exp',
          'hmax',         1000,
          'splitcorners', 1)

        # Quick diagnostic: mesh span in x-direction
        print(md.mesh.x.max() -  md.mesh.x.min())

    elif modelnum == 2 or modelnum == 4:
        md = bamg(model, domain='./Domain.exp', hmax=2000, splitcorners=1)
    elif modelnum == 5 or modelnum == 6:
        md = bamg(model, domain='./Domain.exp', hmax=500,  splitcorners=1)
    elif modelnum == 7 or modelnum == 8:
        md = bamg(model, domain='./Domain.exp', hmax=200,  splitcorners=1)
    else:
        raise RuntimeError("Model not supported yet")

    md.miscellaneous.name = 'MISMIP_' + modelname
    org.savemodel(md)

# ============================================================================
#  STEP2: Parameterization
# ============================================================================
if org.perform('Parameterization'):
    md = org.loadmodel('Mesh_generation')

    # Blank mask (no domain partitioning yet) — customize as needed
    md = setmask(md, '', '')

    # Apply problem-specific parameters (friction, geometry, BCs, etc.)
    # The file Mismip.py should fill in md.* fields accordingly.
    md = parameterize(md, './Mismip.py')

    org.savemodel(md)

# ============================================================================
#  STEP3: Transient_Steadystate (very long relaxation on SSA)
# ============================================================================
if org.perform('Transient_Steadystate'):
    md = org.loadmodel('Parameterization')
    md = setflowequation(md, 'SSA', 'all')

    # For "coulomb" variants, switch friction law and parameters
    if (modelnum == 3) or (modelnum == 4) or (modelnum == 6):
        md.friction = frictioncoulomb()
        md.friction.coefficient = math.sqrt(3.160e6) * np.ones(md.mesh.numberofvertices)
        md.friction.coefficientcoulomb = math.sqrt(0.5) * np.ones(md.mesh.numberofvertices)
        md.friction.p = 3 * np.ones(md.mesh.numberofelements)
        md.friction.q = np.zeros(md.mesh.numberofelements)

    # Multi-kyr relaxation to a quasi steady geometry
    md.timestepping.time_step = 1
    md.timestepping.final_time = 200000
    md.settings.output_frequency = 2000
    md.settings.checkpoint_frequency = 2000

    # Conservative nonlinear solver settings for robust relaxation
    md.stressbalance.maxiter = 30
    md.stressbalance.abstol = float('nan')
    md.stressbalance.restol = 1

    md.verbose = verbose('solution', True, 'module', True, 'convergence', True)
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_Tss1'
    md.settings.solver_residue_threshold = float('nan')

    # Transient (geometry evolves) run; `loadonly` submits or prepares run
    solutiontype = 'tr'
    md = solve(md, solutiontype, 'loadonly', loadonly, 'lock', lock, 'runtimename', False)

    # Optional NetCDF export for portability/inspection
    export_netCDF(md, "./mismip_1km_viscous_Transient_steadystate.nc")

    org.savemodel(md)

# ============================================================================
#  STEP4: Transient_Steadystate_remesh (reuse relaxed state, then remesh)
# ============================================================================
if org.perform('Transient_Steadystate_remesh'):
    md = org.loadmodel('Parameterization')

    # Example: seed from another relaxed solution at a different resolution
    # (adjust path to an existing model in your workspace)
    md2 = loadmodel('Models_2km_viscous/mismip_2km_viscous_Transient_steadystate2.nc')

    # Transfer final solution mesh to current model (simple overwrite approach)
    md2.results.TransientSolution[-1].MeshX = md.mesh.x
    md2.results.TransientSolution[-1].MeshY = md.mesh.y
    md2.results.TransientSolution[-1].MeshElements = md.mesh.elements

    # Remesh guided by your parameter file (e.g., hmax fields, refinement logic)
    md = remesh(md2, './Mismip.py')

    # Reapply flow equation and (if needed) Coulomb friction
    md = setflowequation(md, 'SSA', 'all')
    if (modelnum == 3) or (modelnum == 4) or (modelnum == 6):
        md.friction = frictioncoulomb()
        md.friction.coefficient = math.sqrt(3.160e6) * np.ones(md.mesh.numberofvertices)
        md.friction.coefficientcoulomb = math.sqrt(0.5) * np.ones(md.mesh.numberofvertices)
        md.friction.p = 3 * np.ones(md.mesh.numberofelements)
        md.friction.q = np.zeros(md.mesh.numberofelements)

    # Shorter relaxation/checkpointing for the remeshed state
    md.timestepping.time_step = 1
    md.timestepping.final_time = 200000
    md.settings.output_frequency = 5000
    md.settings.checkpoint_frequency = 5000

    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = float('nan')
    md.stressbalance.restol = 1

    md.verbose = verbose('convergence', False, 'solution', False)
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_Tssr'
    md.settings.solver_residue_threshold = float('nan')
    md.settings.waitonlock = 0

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP5: Transient_steadystate2 (restart from STEP3 end-state)
# ============================================================================
if org.perform('Transient_steadystate2'):
    md = org.loadmodel('Transient_Steadystate')

    md = setflowequation(md, 'SSA', 'all')

    # Re-initialize from the last saved transient state
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    # Shorter follow-on run to tighten steady state
    md.timestepping.time_step = 1
    md.timestepping.final_time = 10000
    md.settings.output_frequency = 1000
    md.settings.checkpoint_frequency = 5000

    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = float('nan')
    md.stressbalance.restol = 1

    md.verbose = verbose('solution', True, 'module', True, 'convergence', True)
    md.cluster = cluster
    md.settings.waitonlock = 1
    md.miscellaneous.name = 'MISMIP_' + modelname + '_Tss2'

    solutiontype = 'tr'
    md = solve(md, solutiontype, 'loadonly', loadonly, 'lock', lock)
    # export_netCDF(...)  # optional

# ============================================================================
#  STEP6: Transient_steadystate3 (another tightening pass)
# ============================================================================
if org.perform('Transient_steadystate3'):
    md = org.loadmodel('Transient_steadystate2')
    md = setflowequation(md, 'SSA', 'all')

    # Reuse the last transient state as initial condition
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    md.timestepping.time_step = 1
    md.timestepping.final_time = 200000
    md.settings.output_frequency = 6000

    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = float('nan')
    md.stressbalance.restol = 1

    md.verbose = verbose('solution', True, 'module', True, 'convergence', True)
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_Tss3'
    md.settings.solver_residue_threshold = float('nan')
    md.settings.waitonlock = 0

    solutiontype = 'tr'
    md = solve(md, solutiontype, 'loadonly', loadonly, 'lock', lock)
    org.savemodel(md)

# ============================================================================
#  STEP7: Transient_steadystate4
# ============================================================================
if org.perform('Transient_steadystate4'):
    md = org.loadmodel('Transient_steadystate3')
    md = setflowequation(md, 'SSA', 'all')

    # Reuse last state
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    md.timestepping.time_step = 1
    md.timestepping.final_time = 200000
    md.settings.output_frequency = 6000

    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = float('nan')
    md.stressbalance.restol = 1

    md.verbose = verbose('convergence', False, 'solution', True)
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_Tss4'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP8: Transient_steadystate5
# ============================================================================
if org.perform('Transient_steadystate5'):
    md = org.loadmodel('Transient_steadystate4')
    md = setflowequation(md, 'SSA', 'all')

    # Reuse last state
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    md.timestepping.time_step = 1
    md.timestepping.final_time = 200000
    md.settings.output_frequency = 6000

    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = float('nan')
    md.stressbalance.restol = 1

    md.verbose = verbose('convergence', False, 'solution', True)
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_Tss5'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP9: Transient_extrude (prepare 3D for HO/FS branches)
# ============================================================================
if org.perform('Transient_extrude'):
    md = org.loadmodel('Transient_steadystate3')

    # Initialize from last transient state
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    # Extrude to 3D: 10 layers, geometric stretching 1.1
    md = md.extrude(10, 1.1)

    # Switch to HO for subsequent runs; keep physics lightweight for now
    md = setflowequation(md, 'HO', 'all')
    md.transient.isthermal = 0
    md.transient.issmb = 0
    md.initialization.temperature[:] = 273.0

    org.savemodel(md)

# ============================================================================
#  STEP10: GlenSSA — Glen's flow law under SSA
# ============================================================================
if org.perform('GlenSSA'):
    # Choose base state: 500m meshes can reuse 2D; others take extruded 3D then collapse
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    md = md.collapse()  # collapse 3D back to 2D
    md = setflowequation(md, 'SSA', 'all')

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_GSSA'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP11: GlenMOLHO — Glen + MOLHO approx with shear BCs
# ============================================================================
if org.perform('GlenMOLHO'):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'VxShear','VyShear','VxBase','VyBase','VxSurface','VySurface',
        'VxAverage','VyAverage','StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    md = md.collapse()
    md = setflowequation(md, 'MOLHO', 'all')

    # Example solver configuration; add block preconditioner for robustness
    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.toolkits = toolkits.addoptions(md.toolkits,'StressbalanceAnalysis',bcgslbjacobioptions())

    # Apply MOLHO shear boundary conditions (helper provided by ISSM)
    md = SetMOLHOBC(md)

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_GMOLHO'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP12: GlenHO — Glen under HO (3D solved but often collapsed for speed)
# ============================================================================
if org.perform('GlenHO'):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md = setflowequation(md, 'HO', 'all')
    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRatexz','StrainRateyz','StrainRatezz',
        'StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # Optional: block CG + Jacobi preconditioner
    # md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_GHO'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP13: GlenFS — short FS sanity run with free surface & GL migration
# ============================================================================
if org.perform('GlenFS'):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_steadystate3')

    # Initialize from a steady state before extruding
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    md = md.extrude(10, 1.1)
    md = setflowequation(md, 'HO', 'all')
    md.transient.isthermal = 0
    md.transient.issmb = 0
    md.initialization.temperature[:] = 273.0

    # Switch to FS for a very short diagnostic transient
    md = setflowequation(md, 'FS', 'all')
    md.stressbalance.shelf_dampening = 1
    md.masstransport.isfreesurface = 1
    md.transient.isgroundingline = 1
    md.groundingline.migration='Contact'  # alternative: 'SoftMigration'

    md.timestepping.time_step = 0.00001
    md.timestepping.final_time = 0.0001
    md.settings.output_frequency = 1

    md.toolkits = toolkits.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())
    md.flowequation.fe_FS = 'TaylorHood'

    md.stressbalance.maxiter = 20
    md.stressbalance.restol = 0.5
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_GFS'

    solutiontype = 'tr'
    md = solve(md, solutiontype, 'loadonly', loadonly, 'lock', lock, 'runtimename', False)
    # org.savemodel(md)  # optional

# ============================================================================
#  STEP14: GlenESSA — SSA with multiplicative enhancement E
# ============================================================================
if org.perform('GlenESSA'):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    md = collapse(md)
    md = setflowequation(md, 'SSA', 'all')

    # Apply Glen enhancement E
    md.materials = matenhancedice(md.materials)
    md.materials.rheology_E = 5.0 * np.ones(md.mesh.numberofvertices)

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_ESSA'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP15: GlenEMOLHO — MOLHO with enhancement E
# ============================================================================
if org.perform('GlenEMOLHO'):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'VxShear','VyShear','VxBase','VyBase','VxSurface','VySurface',
        'VxAverage','VyAverage','StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    md = collapse(md)
    md = setflowequation(md, 'MOLHO', 'all')

    md.materials = matenhancedice(md.materials)
    md.materials.rheology_E[:] = 5.0

    md = SetMOLHOBC(md)

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_GEMOLHO'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP16: GlenEHO — HO with enhancement E
# ============================================================================
if org.perform('GlenEHO'):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md = setflowequation(md, 'HO', 'all')
    md.materials = matenhancedice(md.materials)
    md.materials.rheology_E = 5.0 * np.ones(md.mesh.numberofvertices)

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRatexz','StrainRateyz','StrainRatezz',
        'StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = True
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_GEHO'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP17: GlenEFS — FS with enhancement E (very short diagnostic)
# ============================================================================
if org.perform('GlenEFS'):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md = setflowequation(md, 'FS', 'all')
    md.materials = matenhancedice(md.materials)
    md.materials.rheology_E[:] = 5.0

    md.stressbalance.shelf_dampening = 1
    md.masstransport.isfreesurface = 1
    md.groundingline.migration = 'Contact'

    md.timestepping.time_step = 0.001
    md.timestepping.final_time = 0.002
    md.settings.output_frequency = 1

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = True
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_GEFS'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP (ESTAR variants): ESTARSSA / ESTARMOLHO / ESTARHO
# ============================================================================
# These blocks demonstrate ESTAR (strain-rate–dependent) enhancements.
# NOTE: The original snippet had small typos (e.g., `md,materials`). We keep
# structure but annotate potential corrections in comments.

if org.perform("ESTARSSA"):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md = setflowequation(md, 'SSA', 'all')
    md.materials = matenhancedice(md.materials)  # or matestar(md.materials) if available

    # Two-parameter ESTAR: Es (shear) and Ec (compression) multipliers
    md.materials.rheology_Es = 5.0 * np.ones(md.mesh.numberofvertices)
    # Potential fix: `md.materials`, not `md,materials` (typo in original)
    md.materials.rheology_Ec = (3.0/8.0) * md.materials.rheology_Es

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective','LambdaS','Epsprime',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset']

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_ESTARSSA'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

if org.perform("ESTARMOLHO"):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'VxShear','VyShear','VxBase','VyBase','VxSurface','VySurface','VxAverage','VyAverage',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective','LambdaS','Epsprime',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset']

    md = collapse(md)
    md = setflowequation(md, 'MOLHO', 'all')

    # If your ISSM has a dedicated ESTAR material helper, use that instead:
    # md.materials = matestar(md.materials)
    md.materials = matenhancedice(md.materials)
    md.materials.rheology_Es = 5.0 * np.ones(md.mesh.numberofvertices)
    md.materials.rheology_Ec = (3.0/8.0) * md.materials.rheology_Es

    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    # Optional: tighter solver for very fine meshes
    # (original logic referenced `res`; wire up if you track resolution)
    # if res <= 1000:
    #     md.toolkits = toolkits()
    #     md.toolkits = addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())
    #     md.settings.solver_residue_threshold = np.nan
    #     md.stressbalance.maxiter = 50
    #     md.stressbalance.restol = 1e-4
    #     md.stressbalance.reltol = np.nan
    #     md.stressbalance.abstol = np.nan

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_ESTARMOLHO'

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# STEP18 ---------------------------------------------------------------------
if org.perform("ESTARHO"):
    if modelnum == 5 or modelnum == 6:
        md = org.loadmodel('Transient_Steadystate')
    else:
        md = org.loadmodel('Transient_extrude')

    md = setflowequation(md, 'HO', 'all')

    # As above: use matestar if available; else reuse matenhancedice as a stand-in
    # to carry Es/Ec fields for your postprocessing/diagnostics.
    # md.materials = matestar(md.materials)
    md.materials = matenhancedice(md.materials)

    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRatexz','StrainRateyz','StrainRatezz',
        'StrainRateeffective','LambdaS',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset']

    md.materials.rheology_Es = 5.0 * np.ones(md.mesh.numberofvertices)
    md.materials.rheology_Ec = (3.0/8.0) * md.materials.rheology_Es

    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # md.toolkits = addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())
    md.stressbalance.maxiter = 30
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = float('nan')

    md.verbose.convergence = False
    md.cluster = cluster
    md.miscellaneous.name = 'MISMIP_' + modelname + '_EFS'  # label kept from original

    solutiontype = 'tr'
    md = solve(md, solutiontype)
    org.savemodel(md)

# ============================================================================
#  STEP19: Simple comparative analysis/plotting (placeholder)
# ============================================================================
if org.perform('analyse'):
    # Load different approximations to compare surface/shared fields
    mdgs = org.loadmodel('GlenSSA')
    mdgm = org.loadmodel('GlenMOLHO')
    mdgh = org.loadmodel('GlenHO')
    mdes = org.loadmodel('ESTARSSA')  # or GlenESSA
    mdeh = org.loadmodel('ESTARHO')   # or GlenEHO

    mdgsV = mdgs.results.TransientSolution[-1].Vel
    mdgmV = mdgm.results.TransientSolution[-1].Vel
    mdghV = mdgh.results.TransientSolution[-1].Vel
    mdesV = mdes.results.TransientSolution[-1].Vel
    mdehV = mdeh.results.TransientSolution[-1].Vel

    # NOTE: The original used function-call syntax on arrays; in NumPy use []
    # mdehV[mdeh.mesh.vertexonsurface==1] would subset surface nodes.
    # The call below keeps the original structure for visibility and should be
    # adapted to your plotting utilities.
    plotmodel(
        mdgs,
        'data', mdgsV,
        'data', mdgmV,
        'data', mdghV,  # consider mdghV[mdeh.mesh.vertexonsurface==1]
        'data', mdesV,
        'data', mdehV,  # consider mdehV[mdeh.mesh.vertexonsurface==1]
        'nlines', 5,
        'ncols', 1,
        'caxis#all', [1, 1000]
    )
    plt.show()

print("Done with the MISMIP script in Python.")
