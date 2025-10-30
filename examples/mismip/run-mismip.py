import os
import sys
import argparse
import numpy as np

from loadmodel import loadmodel
from solve import solve
from model import *
from bamg import bamg
from setmask import setmask
from parameterize import parameterize
from export_netCDF import export_netCDF
from setflowequation import setflowequation
## TODO: Use gadi rather than gadilocal when merged into main and available in py-tools
from gadi import gadi
# from gadilocal import gadi
from verbose import verbose
from frictioncoulomb import frictioncoulomb
from MatlabFuncs import *

# Check that ISSM_DIR environment variable is set
issm_dir = os.getenv('ISSM_DIR')
if not issm_dir:
    print("Error: ISSM_DIR environment variable is not set correcttly. Try re-running:\n"
          "   module use /g/data/vk83/prerelease/modules\n"
          "   module load access-issm/pr26-1")
    sys.exit(1)

# Parse command-line arguments
parser = argparse.ArgumentParser(description = 'Run MISMIP experiments using ACCESS-ISSM.')
parser.add_argument('--project_code', type = str, required = True, help = 'Project code for job submission')
parser.add_argument('--steps', type = int, nargs = '+', required = True, help = 'List of steps to run the experiment')
parser.add_argument('--execution_dir', type = str, required = True, help = 'Directory to execute the job in')
parser.add_argument('--storage', type = str, required = True, help = 'PBS storage location for job submission')

parser.add_argument('--model_num', type = int, default = 1, help = 'Model number to use for the experiment')
parser.add_argument('--load_only', type = bool, default = False, help = 'Should results be loaded only without model execution?')
parser.add_argument('--num_nodes', type = int, default = 1, help = 'Number of nodes to use for the PBS job')
parser.add_argument('--cpu_node', type = int, default = 32, help = 'Number of CPUs per node to use for the PBS job')
parser.add_argument('--walltime', type = str, default = 2880, help = 'Walltime for the PBS job')
parser.add_argument('--queue', type = str, default = 'normal', help = 'Queue for the PBS job')
## TODO: Update default module load and use for official release
parser.add_argument('--module_use', type = str, nargs = '+', default = ['/g/data/vk83/prerelease/modules'], help = 'Module locations for the PBS job')
parser.add_argument('--module_load', type = str, nargs = '+', default = ['access-issm/pr26-18'], help = 'Modules to load for the PBS job')
parser.add_argument('--memory', type = int, default = 128, help = 'Memory (in GB) to allocate for the PBS job')

args = parser.parse_args()

# Error checks for arguments
if args.model_num < 1 or args.model_num > 8:
    raise RuntimeError("Invalid model number: model_num must be between 1 and 8")
if not all(step in [0, 1, 2, 3] for step in args.steps):
    raise RuntimeError("Invalid steps: steps must be a list containing any of [0, 1, 2, 3]")
if len(args.module_load) != len(args.module_use):
    raise RuntimeError("module_load and module_use must have the same number of entries")


# Set variables from parsed arguments
model_num = args.model_num
model_name = f'mismip_model_{model_num}'
out_dir = os.path.join(args.execution_dir, model_name)
storage_flag = args.storage + '+gdata/vk83'  # Append gdata/vk83 for ISSM executable access
login = os.getlogin()

# Create output directory if it doesn't exist
os.makedirs(out_dir, exist_ok=True)

# Setup cluster configuration
cluster = gadi('name', oshostname(),
               'numnodes', args.num_nodes,
               'cpuspernode', args.cpu_node,
               'memory', args.memory,
               'moduleuse', args.module_use,
               'moduleload', args.module_load,
               'storage', storage_flag,
               'login', login,
               'srcpath', issm_dir,
               'codepath', issm_dir + '/bin',
               'executionpath', out_dir,
               'project', args.project_code,
               'queue', args.queue,
               'time', args.walltime
               )

# --------------------------------------------------------------
# STEP 0: Print configuration settings
# --------------------------------------------------------------
if 0 in args.steps:
    print("=============================================================\n"
            " ACCESS-ISSM MISMIP+ CONFIGURATION SETTINGS \n"
            "=============================================================\n"
            f" ISSM directory:              {issm_dir}\n"
            f" Python executable:           {sys.executable}\n"
            "---------------------------------------------\n"
            f" Project code:                {args.project_code}\n"
            f" User login:                  {login}\n"
            f" Storage locations:           {storage_flag}\n"
            f" Model number:                {args.model_num}\n"
            f" Execution directory:         {args.execution_dir}\n"
            f" Model name:                  {model_name}\n"
            f" Output directory:            {out_dir}\n"
            f" Steps to run:                {args.steps}\n"
            f" Load only:                   {args.load_only}\n"
            f" Walltime:                    {args.walltime}\n"
            f" Queue:                       {args.queue}\n"
            f" Module use locations:        {args.module_use}\n"
            f" Modules to load:             {args.module_load}\n"
            "=============================================================\n")

# --------------------------------------------------------------
# STEP 1: Generate model mesh
# --------------------------------------------------------------
if 1 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 1: Generating mesh for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'mesh'

    # Initialize an empty model
    model = model()

    # Choose hmax based on model_num (different resolutions)
    if model_num == 1 or model_num == 3:
        hmax = 1000
    elif model_num == 2 or model_num == 4:
        hmax = 2000
    elif model_num == 5 or model_num == 6:
        hmax = 500
    elif model_num == 7 or model_num == 8:
        hmax = 200
    else:
        raise RuntimeError("Invalid model number: model_num must be between 1 and 8")

    # Generate the mesh using BAMG
    print(f"Generating mesh for model number {model_num}...")
    md = bamg(model, 'domain', './Domain.exp', 'hmax', hmax, 'splitcorners', 1)

    # Print mesh summary
    print(f"Mesh generated with:\n"
            f"    {md.mesh.numberofvertices} vertices\n"
            f"    {md.mesh.numberofelements} elements")

    # Set model name
    md.miscellaneous.name = model_name

    # Save the model
    print(f"Saving mesh to {out_dir}...")
    export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 2: Parameterise the model
# --------------------------------------------------------------
if 2 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 2: Parameterising model for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'parameterise'

    # Load the mesh created in Step 1
    print(f"Loading mesh from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_mesh.nc'))
        
    # Set mask (assume all ice is grounded)
    print(f"Parameterising model for model number {model_num}...")
    md = setmask(md, '', '')

    # Set model parameters
    md = parameterize(md, './mismip_param.py')

    # Save the parameterised model
    print(f"Saving parameterised model to {out_dir}...")
    export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 3: Run Transient Steady-State Simulation
# --------------------------------------------------------------
if 3 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 3: Running Transient Steady-State Simulation for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_steady_state'

    # Load the parameterised model from Step 2
    print(f"Loading parameterised model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_parameterise.nc'))

    # Set flow equation to SSA
    print(f"Setting flow equation to SSA for model number {model_num}...")
    md = setflowequation(md, 'SSA', 'all')

    # For "coulomb" variants, switch friction law and parameters
    if (model_num == 3) or (model_num == 4) or (model_num == 6):
        print(f"Adjusting friction law to Coulomb for model number {model_num}...")
        md.friction = frictioncoulomb()
        md.friction.coefficient = math.sqrt(3.160e6) * np.ones(md.mesh.numberofvertices)
        md.friction.coefficientcoulomb = math.sqrt(0.5) * np.ones(md.mesh.numberofvertices)
        md.friction.p = 3 * np.ones(md.mesh.numberofelements)
        md.friction.q = np.zeros(md.mesh.numberofelements)

    # Set transient simulation parameters
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1
    md.timestepping.final_time = 20 # TODO: Set to 200000
    md.settings.output_frequency = 5 # TODO: Set to 2000
    md.settings.checkpoint_frequency = 5 # TODO: Set to 2000

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 30
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.verbose = verbose('solution', True, 'module', True, 'convergence', True)
    md.cluster = cluster
    md.miscellaneous.name = model_name + '_transient_steady_state'
    md.settings.solver_residue_threshold = np.nan
    md.settings.waitonlock = 0

    # Solve/load the transient steady-state simulation
    if args.load_only:
        print(f"Loading transient steady-state results from {out_dir}...")
    else:
        print(f"Solving transient steady-state simulation for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving transient steady-state results to {out_dir}...")
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 4: Run Transient Steady-State Remeshing Simulation
# --------------------------------------------------------------
if 4 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 4: Running Transient Steady-State Remeshing Simulation for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_steady_state_remeshing'

    # Load the parameterised model from Step 2
    print(f"Loading parameterised model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_parameterise.nc'))

    ## TODO: How do we handle this optional remeshing model step? It requires a model path to be specified


if 5 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 5: Running second Transient Steady-State Simulation for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_steady_state_2'

    # Load the transient steady-state results from Step 3
    print(f"Loading transient steady-state results from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_steady_state.nc'))

    # Set flow equation to SSA
    # NOTE: Only necessary if a remeshed model is loaded, otherwise retained from Step 3
    md = setflowequation(md, 'SSA', 'all')

    # Re-initialize from the last saved transient state
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    # Set transient simulation parameters (shorter simulation)
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1
    md.timestepping.final_time = 10 ## TODO: Set to 10000
    md.settings.output_frequency = 2 ## TODO: Set to 1000
    md.settings.checkpoint_frequency = 2 ## TODO: Set to 5000

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.verbose = verbose('solution', True, 'module', True, 'convergence', True)
    md.cluster = cluster
    md.settings.waitonlock = 1

    # Solve/load the second transient steady-state simulation
    if args.load_only:
        print(f"Loading second transient steady-state results from {out_dir}...")
    else:
        print(f"Solving second transient steady-state simulation for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving second transient steady-state results to {out_dir}...")
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))
