import os
import sys
import argparse
import numpy as np

temp = False

if temp:
    os.environ['ISSM_DIR'] = '/g/data/vk83/apps/spack/0.22/release/linux-rocky8-x86_64/gcc-13.2.0/issm-git.2025.11.24_2025.11.24-pd55xlx56v5vuno2lshenunfddfvupnr'
    sys.path.append(os.getenv('ISSM_DIR') + '/lib')
    sys.path.append(os.getenv('ISSM_DIR') + '/python-tools.zip')
    os.chdir('/scratch/tm70/lb9857/ACCESS-ISSM/examples/mismip')

from loadmodel import loadmodel
from solve import solve
from model import *
from bamg import bamg
from setmask import setmask
from parameterize import parameterize
from export_netCDF import export_netCDF
from setflowequation import setflowequation
from gadi import gadi
from verbose import verbose
from frictioncoulomb import frictioncoulomb
from MatlabFuncs import *

# Check that ISSM_DIR environment variable is set
issm_dir = os.getenv('ISSM_DIR')
if not issm_dir:
    print("Error: ISSM_DIR environment variable is not set correcttly. Try re-running:\n"
          "   module use /g/data/vk83/modules\n"
          "   module load access-issm/2025.11.0")
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
parser.add_argument('--module_use', type = str, nargs = '+', default = ['/g/data/vk83/modules'], help = 'Module locations for the PBS job')
parser.add_argument('--module_load', type = str, nargs = '+', default = ['access-issm/2025.11.0'], help = 'Modules to load for the PBS job')
parser.add_argument('--memory', type = int, default = 128, help = 'Memory (in GB) to allocate for the PBS job')

args = parser.parse_args()

if temp:
    args = argparse.Namespace(
        project_code='tm70',
        steps=[5],
        execution_dir='/scratch/tm70/lb9857/ACCESS-ISSM/examples/mismip/execution',
        storage='scratch/tm70',
        model_num=1,
        load_only=False,
        num_nodes=1,
        cpu_node=32,
        walltime=2880,
        queue='normal',
        module_use=['/g/data/vk83/modules'],
        module_load=['access-issm/2025.11.0'],
        memory=128
    )

# Error checks for arguments
if args.model_num < 1 or args.model_num > 8:
    raise RuntimeError("Invalid model number: model_num must be between 1 and 8")
if not all(step in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] for step in args.steps):
    raise RuntimeError("Invalid steps: steps must be a list containing any of [0, 1, 2, 3, 4, 5, 6, 7, 8]")
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
    md.miscellaneous.name = f'{model_name}_{step_name}'

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
    md.miscellaneous.name = f'{model_name}_{step_name}'

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
    md.timestepping.final_time = 200000 # TODO: Set to 200000
    md.settings.output_frequency = 2000 # TODO: Set to 2000
    md.settings.checkpoint_frequency = 2000 # TODO: Set to 2000 ## DO WE NEED THIS? IF IT DOESN'T EXCEED WALL TIME, WE DON'T NEED THE RESTARTS

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 30
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
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
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 4: Run second Transient Steady-State Simulation
# --------------------------------------------------------------
if 4 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 4: Running second Transient Steady-State Simulation for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_steady_state_2'

    # Load the transient steady-state results from Step 3
    print(f"Loading transient steady-state results from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_steady_state.nc'))

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
    md.settings.checkpoint_frequency = 2 ## TODO: Set to 5000 ## DO WE NEED THIS? IF IT DOESN'T EXCEED WALL TIME, WE DON'T NEED THE RESTARTS

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the second transient steady-state simulation
    if args.load_only:
        print(f"Loading second transient steady-state results from {out_dir}...")
    else:
        print(f"Solving second transient steady-state simulation for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving second transient steady-state results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 5: Run third Transient Steady-State Simulation
# --------------------------------------------------------------
if 5 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 5: Running third Transient Steady-State Simulation for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_steady_state_3'

    # Load the transient steady-state results from Step 4
    print(f"Loading transient steady-state results from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_steady_state_2.nc'))

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
    md.timestepping.final_time = 10 ## TODO: Set to 200000
    md.settings.output_frequency = 2 ## TODO: Set to 6000
    md.settings.checkpoint_frequency = 2 ## TODO: Set to 6000 ## DO WE NEED THIS? IF IT DOESN'T EXCEED WALL TIME, WE DON'T NEED THE RESTARTS

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the third transient steady-state simulation
    if args.load_only:
        print(f"Loading third transient steady-state results from {out_dir}...")
    else:
        print(f"Solving third transient steady-state simulation for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving third transient steady-state results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 6: Run fourth Transient Steady-State Simulation
# --------------------------------------------------------------
if 6 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 6: Running fourth Transient Steady-State Simulation for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_steady_state_4'

    # Load the transient steady-state results from Step 5
    print(f"Loading transient steady-state results from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_steady_state_3.nc'))

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
    md.timestepping.final_time = 10 ## TODO: Set to 200000
    md.settings.output_frequency = 2 ## TODO: Set to 6000
    md.settings.checkpoint_frequency = 2 ## TODO: Set to 6000 ## DO WE NEED THIS? IF IT DOESN'T EXCEED WALL TIME, WE DON'T NEED THE RESTARTS

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the fourth transient steady-state simulation
    if args.load_only:
        print(f"Loading fourth transient steady-state results from {out_dir}...")
    else:
        print(f"Solving fourth transient steady-state simulation for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving fourth transient steady-state results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 7: Run fifth Transient Steady-State Simulation
# --------------------------------------------------------------
if 7 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 7: Running fifth Transient Steady-State Simulation for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_steady_state_5'

    # Load the transient steady-state results from Step 6
    print(f"Loading transient steady-state results from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_steady_state_4.nc'))

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
    md.timestepping.final_time = 10 ## TODO: Set to 200000
    md.settings.output_frequency = 2 ## TODO: Set to 6000
    md.settings.checkpoint_frequency = 2 ## TODO: Set to 6000 ## DO WE NEED THIS? IF IT DOESN'T EXCEED WALL TIME, WE DON'T NEED THE RESTARTS

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the fifth transient steady-state simulation
    if args.load_only:
        print(f"Loading fifth transient steady-state results from {out_dir}...")
    else:
        print(f"Solving fifth transient steady-state simulation for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving fifth transient steady-state results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 8: Extrude mesh to 3D
# --------------------------------------------------------------
if 8 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 8: Extruding mesh to 3D for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'transient_extude'

    # Load the transient steady-state results from Step 5
    print(f"Loading transient steady-state (5) results from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_steady_state_5.nc'))

    # Re-initialize from the last saved transient state
    md.initialization.vx = md.results.TransientSolution[-1].Vx
    md.initialization.vy = md.results.TransientSolution[-1].Vy
    md.initialization.vel = md.results.TransientSolution[-1].Vel
    md.geometry.thickness = md.results.TransientSolution[-1].Thickness
    md.geometry.base = md.results.TransientSolution[-1].Base
    md.geometry.surface = md.results.TransientSolution[-1].Surface
    md.mask.ocean_levelset = md.results.TransientSolution[-1].MaskOceanLevelset

    # Extrude to 3D: 10 layers, geometric stretching 1.1
    md = md.extrude(10, 1.1)

    # Switch to HO for subsequent runs
    md = setflowequation(md, 'HO', 'all')
    md.transient.isthermal = 0
    md.transient.issmb = 0
    md.initialization.temperature[:] = 273.0

    # Save model
    print(f"Saving extruded 3D model to {out_dir}...")
    export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 9: EXPERIMENT - Glen_SSA
# --------------------------------------------------------------
if 9 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 9: Glen_SSA experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'Glen_SSA'

    # Load the extruded model from Step 8
    print(f"Loading the 'transient_extrude' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_extrude.nc'))

    # Define requested outputs
    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    # Collapse model back to 2D and set SSA flow equation
    print(f"Collapsing model to 2D and setting SSA flow equation for model number {model_num}...")
    md = md.collapse()
    md = setflowequation(md, 'SSA', 'all')

    # Set transient simulation parameters (shorter simulation)
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_SSA experiment
    if args.load_only:
        print(f"Loading Glen_SSA experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_SSA experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_SSA experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 10: EXPERIMENT - Glen_HO
# --------------------------------------------------------------
if 10 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 10: Glen_HO experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'Glen_HO'

    # Load the extruded model from Step 8
    print(f"Loading the 'transient_extrude' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_extrude.nc'))

    # Define requested outputs
    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRatexz','StrainRateyz','StrainRatezz',
        'StrainRateeffective', 'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    # Set HO flow equation
    print(f"Setting HO flow equation for model number {model_num}...")
    md = setflowequation(md, 'HO', 'all')

    # Adjust toolkit options for better convergence
    print(f"Adjusting the toolkits options for model number {model_num}...")
    md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    # Set transient simulation parameters
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_HO experiment
    if args.load_only:
        print(f"Loading Glen_HO experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_HO experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_HO experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 11: EXPERIMENT - Glen_FS
# --------------------------------------------------------------
if 11 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 11: Glen_FS experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    raise NotImplementedError("Glen_FS experiment has not been tested yet.\n"
                              "Initial implementation is provided in `mismip/run-mismip.py`, but community testing is required to verify correctness.\n"
                              "If you conduct testing, please contribute your results back to ACCESS-NRI.\n")

    step_name = 'Glen_FS'

    # Load the output from Glen_HO experiment (Step 10) for improved initial conditions
    print(f"Loading the 'Glen_HO' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_Glen_HO.nc'))

    # Define requested outputs
    md.transient.requested_outputs = ['default']

    # Set FS flow equation
    print(f"Setting FS flow equation for model number {model_num}...")
    md = setflowequation(md, 'FS', 'all')
    md.flowequation.fe_FS = 'TaylorHood'

    # Adjust toolkit options for better convergence
    print(f"Adjusting the toolkits options for model number {model_num}...")
    md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    # Set transient simulation parameters (conduct a short test run)
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 0.00001
    md.timestepping.final_time = 0.0001
    md.settings.output_frequency = 1
    md.transient.isgroundingline = 1
    md.transient.ismasstransport = 1
    md.groundingline.migration = 'Contact'

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 20
    md.stressbalance.restol = 0.5
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan
    md.stressbalance.shelf_dampening = 1
    md.masstransport.isfreesurface = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_FS experiment
    if args.load_only:
        print(f"Loading Glen_FS experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_FS experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_FS experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 12: EXPERIMENT - Glen_E_SSA
# --------------------------------------------------------------
if 12 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 12: Glen_E_SSA experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'Glen_E_SSA'

    # Load the extruded model from Step 8
    print(f"Loading the 'transient_extrude' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_extrude.nc'))

    # Define requested outputs
    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRateeffective',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    # Collapse model back to 2D and set SSA flow equation
    print(f"Collapsing model to 2D and setting SSA flow equation for model number {model_num}...")
    md = md.collapse()
    md = setflowequation(md, 'SSA', 'all')

    # Apply Glen E enhancement
    print(f"Applying Glen E enhancement for model number {model_num}...")
    md.materials = matenhancedice(md.materials)
    md.materials.rheology_E[:] = 5.0

    # Set transient simulation parameters (shorter simulation)
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_E_SSA experiment
    if args.load_only:
        print(f"Loading Glen_E_SSA experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_E_SSA experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_E_SSA experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 13: EXPERIMENT - Glen_E_HO
# --------------------------------------------------------------
if 13 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 13: Glen_E_HO experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'Glen_E_HO'

    # Load the extruded model from Step 8
    print(f"Loading the 'transient_extrude' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_extrude.nc'))

    # Define requested outputs
    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRatexz','StrainRateyz','StrainRatezz',
        'StrainRateeffective', 'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    # Set HO flow equation
    print(f"Setting HO flow equation for model number {model_num}...")
    md = setflowequation(md, 'HO', 'all')

    # Apply Glen E enhancement
    print(f"Applying Glen E enhancement for model number {model_num}...")
    md.materials = matenhancedice(md.materials)
    md.materials.rheology_E[:] = 5.0

    # Adjust toolkit options for better convergence
    print(f"Adjusting the toolkits options for model number {model_num}...")
    md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    # Set transient simulation parameters
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_E_HO experiment
    if args.load_only:
        print(f"Loading Glen_E_HO experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_E_HO experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_E_HO experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 14: EXPERIMENT - Glen_E_FS
# --------------------------------------------------------------
if 14 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 14: Glen_E_FS experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    raise NotImplementedError("Glen_E_FS experiment has not been tested yet.\n"
                              "Initial implementation is provided in `mismip/run-mismip.py`, but community testing is required to verify correctness.\n"
                              "If you conduct testing, please contribute your results back to ACCESS-NRI.\n")

    step_name = 'Glen_E_FS'

    # Load the output from Glen_E_HO experiment (Step 13) for improved initial conditions
    print(f"Loading the 'Glen_E_HO' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_Glen_E_HO.nc'))

    # Define requested outputs
    md.transient.requested_outputs = ['default']

    # Set FS flow equation
    print(f"Setting FS flow equation for model number {model_num}...")
    md = setflowequation(md, 'FS', 'all')
    md.flowequation.fe_FS = 'TaylorHood'

    # Apply Glen E enhancement
    print(f"Applying Glen E enhancement for model number {model_num}...")
    md.materials = matenhancedice(md.materials)
    md.materials.rheology_E[:] = 5.0

    # Adjust toolkit options for better convergence
    print(f"Adjusting the toolkits options for model number {model_num}...")
    md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    # Set transient simulation parameters (conduct a short test run)
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 0.00001
    md.timestepping.final_time = 0.0001
    md.settings.output_frequency = 1
    md.transient.isgroundingline = 1
    md.transient.ismasstransport = 1
    md.groundingline.migration = 'Contact'

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 20
    md.stressbalance.restol = 0.5
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan
    md.stressbalance.shelf_dampening = 1
    md.masstransport.isfreesurface = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_E_FS experiment
    if args.load_only:
        print(f"Loading Glen_E_FS experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_E_FS experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_E_FS experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 15: EXPERIMENT - Estar_SSA
# --------------------------------------------------------------
if 15 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 15: Estar_SSA experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'Estar_SSA'

    # Load the extruded model from Step 8
    print(f"Loading the 'transient_extrude' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_extrude.nc'))

    # Define requested outputs
    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea', 'StrainRatexx',
        'StrainRatexy','StrainRateyy','StrainRateeffective', 'LambdaS','Epsprime',
        'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    # Collapse model back to 2D and set SSA flow equation
    print(f"Collapsing model to 2D and setting SSA flow equation for model number {model_num}...")
    md = md.collapse()
    md = setflowequation(md, 'SSA', 'all')

    # Apply Estar flow law
    print(f"Applying Estar flow law for model number {model_num}...")
    md.materials = matestar(md.materials)
    md.materials.rheology_Es = 5.0 * np.ones(md.mesh.numberofvertices, 1)
    md.materials.rheology_Ec= 3/8 * md.materials.rheology_Es

    # Set transient simulation parameters (shorter simulation)
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Estar_SSA experiment
    if args.load_only:
        print(f"Loading Estar_SSA experiment results from {out_dir}...")
    else:
        print(f"Solving Estar_SSA experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Estar_SSA experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 16: EXPERIMENT - Estar_HO
# --------------------------------------------------------------
if 16 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 16: Estar_HO experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    step_name = 'Estar_HO'

    # Load the extruded model from Step 8
    print(f"Loading the 'transient_extrude' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_transient_extrude.nc'))

    # Define requested outputs
    md.transient.requested_outputs = [
        'default','IceVolume','IceVolumeAboveFloatation','GroundedArea',
        'StrainRatexx','StrainRatexy','StrainRateyy','StrainRatexz','StrainRateyz','StrainRatezz',
        'StrainRateeffective', 'MaskOceanLevelset','IceMaskNodeActivation','MaskIceLevelset'
    ]

    # Set HO flow equation
    print(f"Setting HO flow equation for model number {model_num}...")
    md = setflowequation(md, 'HO', 'all')

    # Apply Estar flow law
    print(f"Applying Estar flow law for model number {model_num}...")
    md.materials = matestar(md.materials)
    md.materials.rheology_Es = 5.0 * np.ones(md.mesh.numberofvertices, 1)
    md.materials.rheology_Ec = 3/8 * md.materials.rheology_Es

    # Adjust toolkit options for better convergence
    print(f"Adjusting the toolkits options for model number {model_num}...")
    md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    # Set transient simulation parameters
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 1.0/12.0
    md.timestepping.final_time = 1000
    md.settings.output_frequency = 600

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 10
    md.stressbalance.restol = 1
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_E_HO experiment
    if args.load_only:
        print(f"Loading Glen_E_HO experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_E_HO experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_E_HO experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))

# --------------------------------------------------------------
# STEP 17: EXPERIMENT - Estar_FS
# --------------------------------------------------------------
if 17 in args.steps:
    print("-------------------------------------------------------------")
    print(f" STEP 17: Estar_FS experiment for MISMIP+ Model {model_num}")
    print("-------------------------------------------------------------")

    raise NotImplementedError("Estar_FS experiment has not been tested yet.\n"
                              "Initial implementation is provided in `mismip/run-mismip.py`, but community testing is required to verify correctness.\n"
                              "If you conduct testing, please contribute your results back to ACCESS-NRI.\n")

    step_name = 'Estar_FS'

    # Load the output from Glen_E_HO experiment (Step 13) for improved initial conditions
    print(f"Loading the 'Estar_HO' model from {out_dir}...")
    md = loadmodel(os.path.join(out_dir, f'{model_name}_Estar_HO.nc'))

    # Define requested outputs
    md.transient.requested_outputs = ['default']

    # Set FS flow equation
    print(f"Setting FS flow equation for model number {model_num}...")
    md = setflowequation(md, 'FS', 'all')
    md.flowequation.fe_FS = 'TaylorHood'

    # Apply Estar flow law
    print(f"Applying Estar flow law for model number {model_num}...")
    md.materials = matestar(md.materials)
    md.materials.rheology_Es = 5.0 * np.ones(md.mesh.numberofvertices, 1)
    md.materials.rheology_Ec = 3/8 * md.materials.rheology_Es

    # Adjust toolkit options for better convergence
    print(f"Adjusting the toolkits options for model number {model_num}...")
    md.toolkits = md.addoptions(md.toolkits, 'StressbalanceAnalysis', bcgslbjacobioptions())

    # Set transient simulation parameters (conduct a short test run)
    print(f"Setting transient simulation parameters for model number {model_num}...")
    md.timestepping.time_step = 0.00001
    md.timestepping.final_time = 0.0001
    md.settings.output_frequency = 1
    md.transient.isgroundingline = 1
    md.transient.ismasstransport = 1
    md.groundingline.migration = 'Contact'

    # Set stress balance parameters
    print(f"Setting stress balance parameters for model number {model_num}...")
    md.stressbalance.maxiter = 20
    md.stressbalance.restol = 0.5
    md.stressbalance.reltol = 0.001
    md.stressbalance.abstol = np.nan
    md.stressbalance.shelf_dampening = 1
    md.masstransport.isfreesurface = 1

    # Set miscellaneous parameters
    print(f"Setting miscellaneous parameters for model number {model_num}...")
    md.cluster = cluster
    md.miscellaneous.name = f'{model_name}_{step_name}'
    md.settings.waitonlock = 0

    # Solve/load the Glen_E_FS experiment
    if args.load_only:
        print(f"Loading Glen_E_FS experiment results from {out_dir}...")
    else:
        print(f"Solving Glen_E_FS experiment for model number {model_num}...")
    md = solve(md, 'transient', 'loadonly', args.load_only, 'runtimename', False)

    # Save the results if loading only
    if args.load_only:
        print(f"Saving Glen_E_FS experiment results to {out_dir}...")
        md.cluster = generic() # TEMPORARY FIX: Clear cluster info before saving to avoid erorr on load due to overzealous login check
        export_netCDF(md, os.path.join(out_dir, f'{model_name}_{step_name}.nc'))