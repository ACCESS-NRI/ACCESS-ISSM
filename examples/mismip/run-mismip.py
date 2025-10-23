import os
import sys
import argparse
import numpy as np

from solve import solve
from model import *
from bamg import bamg
from setmask import setmask
from parameterize import parameterize
from loadmodel import loadmodel
from export_netCDF import export_netCDF
from setflowequation import setflowequation
from gadi import gadi
from verbose import verbose
from frictioncoulomb import frictioncoulomb

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
parser.add_argument('--model_num', type = int, default = 1, help = 'Model number to use for the experiment')
parser.add_argument('--steps', type = int, nargs = '+', required = True, help = 'List of steps to run the experiment')
parser.add_argument('--storage', type = str, required = True, help = 'PBS storage location for job submission')
parser.add_argument('--load_only', type = bool, default = False, help = 'Should results be loaded only without model execution?')
parser.add_argument('--walltime', type = str, default = 2880, help = 'Walltime for the job submission')
parser.add_argument('--queue', type = str, default = 'normal', help = 'Queue to submit the job to')
parser.add_argument('--execution_dir', type = str, default = './execution', help = 'Directory to execute the job in')
parser.add_argument('--num_nodes', type = int, default = 1, help = 'Number of nodes to use for the job')
parser.add_argument('--cpu_node', type = int, default = 32, help = 'Number of CPUs per node to use for the job')

# parser.add_argument('--examples_dir', type = str, default = '~/ACCESS_ISSM/examples', help = 'Source directory for MISMIP+ files')
args = parser.parse_args()

# Set variables from parsed arguments
model_num = args.model_num
model_name = f'mismip_model_{model_num}'
out_dir = os.path.join(args.execution_dir, model_name)
storage_flag = args.storage + '+gdata/vk83'  # Append gdata/vk83 for ISSM executable access
login = os.getlogin()

# Create output directory if it doesn't exist
os.makedirs(out_dir, exist_ok=True)

# Setup cluster configuration
# TODO: Update gadi cluster to operate correctly with all required flags
cluster = gadi('name', 'gadi.nci.org.au',
               'numnodes', args.num_nodes,
               'storage', storage_flag,
               'login', login,
               'srcpath', issm_dir,
               'codepath', issm_dir + '/bin',
               'executionpath', out_dir,
               'project', args.project_code,
               'queue', args.queue,
               'time', args.walltime,
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
    md.timestepping.final_time = 5 # TODO: Set to 200000
    md.settings.output_frequency = 1 # TODO: Set to 2000
    md.settings.checkpoint_frequency = 1 # TODO: Set to 2000

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 30
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.verbose = verbose('solution', True, 'module', True, 'convergence', True)
    md.cluster = cluster
    md.miscellaneous.name = model_name + '_Tss1'
    md.settings.solver_residue_threshold = np.nan

    # Solve the transient steady-state simulation
    print(f"Solving transient steady-state simulation for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving transient steady-state results to {out_dir}...")
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))