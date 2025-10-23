## 2 Prerequisites

!!! warning
    To run {{ model }}, you need to be a member of a project with allocated _Service Units (SU)_. For more information, check [how to join relevant NCI projects](/getting_started/set_up_nci_account#join-relevant-nci-projects).

- **NCI Account**<br>
    Before running {{ model }}, you need to [Set Up your NCI Account](/getting_started/set_up_nci_account).

- **Join NCI projects**<br>
    Join the following projects by requesting membership on their respective NCI project pages:

    - [vk83](https://my.nci.org.au/mancini/project/vk83/join)
    - [xp65](https://my.nci.org.au/mancini/project/xp65/join)

  For more information on joining specific NCI projects, refer to [How to connect to a project](https://opus.nci.org.au/display/Help/How+to+connect+to+a+project).

## 3 Getting started

### 3.1 Obtaining files
An example execution script for the MISMIP+ configuration is available at **[URL]**. To copy these files to your `$HOME` directory, run:

<!-- TODO: Update branch to main -->
```bash
cd
git clone --branch lb/mismip-example https://github.com/ACCESS-NRI/ACCESS-ISSM.git
```

This will create a directory named `ACCESS-ISSM` in your `$HOME` directory. The example execution script is located at `~/ACCESS-ISSM/examples/mismip/run-mismip.py`. This documentation assumes that you do not move these files.

!!! warning
    The execution script uses relative paths. Do not change the file structure of the downloaded files. This documentation assumes that you are _inside_ the `~/ACCESS-ISSM/examples/mismip/` directory when executing all commands.

### 3.1 Workflow overview

Running {{ model }} on _Gadi_ is generally separated into two key components:

* **Model parameterisation** - This step relates to *building* an ISSM model, including defining a model domain and mesh, setting initial conditions, and specifying the model run configuration.
* **Model excution steps** - Subsequent steps relate to various *execution* steps that typically incrementally build upon each previous step.

_Model parameterisation_ is generally less computationally expensive than model execution steps and can tyically be completed simply using a login node on _Gadi_. Depending on the size and configuration of the model, _Model execution_ typically requires a [PBS job](https://opus.nci.org.au/spaces/Help/pages/90308818/4.+PBS+Jobs) submission due to increased computational need.

The general workflow to run {{ model }} is typically separated into individual "steps", contained within a single execution script. The script itself contains all of the model configuration information and handles all model IO. After each step, the existing model state is saved to file as a NetCDF. Below, we provide a brief overview of the key "steps" contained in the example MISMIP+ configuration provided here:

<!-- TODO: Add overview of all execution steps -->
* **Step 1: Model parameterisation #1 (Model domain and mesh configuration)** - 
* **Step 2: Model parameterisation #2 (Initial conditions and configuration)** -
* ...

### 3.2 Setup environment requirements

Interacting with {{ model }} requires the `$ISSM_DIR` and `$PYTHONPATH` environment variables be set to use an appropriate executable and set of Python tools. This is handled automatically when loading the {{ model }} module on _Gadi_. To set these variables in preparation for executing the example script below, run:

<!-- TODO: Update to use official release, not PR26-1 -->
```bash
module use /g/data/vk83/prerelease/modules
module load access-issm/pr26-1
```

In addition, we make use of Python tools to control the modelling workflow and interact with the model files. To prevent the need for all users to maintain individual Python environments, we can leverage the `conda/analysis3` environment maintained by ACCESS-NRI. To load the Python environment, run:

```bash
module use /g/data/xp65/public/modules
module load conda/analysis3
```

### 3.3 Interacting with the execution script
The model execution process is centered around the execution script. Here, we provide a pre-configured execution script for ease of use. Below, we provide key information about interacting with the execution script.

#### 3.3.1 Script arguments
To run {{ model }}, key arguments must be passed to the execution script. Below provides an overview of all available arguments (required and optional) that can be passed to the execution script.

**Required arguments**

The following arguments are required to execute any step(s):

* `--project_code`: No default project code is set. Users must provide a valid _Gadi_ project code for a project they are a member of and that has available compute resources. For example: `--project_code <PROJECT_CODE>`
* `--steps`: No default steps are set. Users must provide one or more steps to execute. Multiple steps should be separated by a space. For example: `--steps 1 2` would run steps 1 and 2. 
* `--execution_dir`: No default execuation directory is set. Users must provide a directory location to store model files. We caution the use of your `$HOME` directory due to storage limitaion. For example: `--execution_dir <EXECUTION_DIR>`
* `--storage`: No default storage flags are set. Users must provide a suitable PBS storage flag to all user required _Gadi_ storage locations (e.g. for your execution directory). For example: `scratch/<PROJECT_CODE>`. Note that `gdata/vk83` is automatically appended, allowing the {{ model }} executable to be loaded within PBS jobs.

**Optional**

The following arguments are optional and primarily adjust default PBS settings set within the execution script:

* `--model_num`: Model configuration to execute (see below). Default = 1.
* `--load_only`: Load results from the specfied step only. No model execution. Default = False.
* `--num_nodes`: Number of nodes to use on _Gadi_. Default = 1.
* `--cpu_node`: Number of CPUs to use per node on _Gadi_. Default = 32.
* `--wall_time`: PBS job walltime (in minutes). Default = 2880.
* `--queue`: Name of _Gadi_ queue to submit the job to. Default = normal.

#### 3.3.2 Example execution submission
The below call provides an example of how the above arguments should be passed to the execution script from your terminal window.

```bash
cd ~/ACCESS-ISSM/examples/mismip/
python run-mismip.py --project_code <PROJECT_CODE --steps 1 2 --execution_dir <EXECUTION_DIR> --load_only False
```

### 3.4 MISMIP Model Configurations
<!-- TODO: Confirm the configuration definitions are correct -->
The provided execution script supports eight different model configurations of differing resolution and friction laws from MISIP+. These configurations are summarised below. The execution script runs Model 1 by default.

| Model number   | Resolution | Friction law| Flow equation |
| -------------- | ---------- | ----------- | ------------- |
| 1              | XX         | XX          | XX            |
| 2              | XX         | XX          | XX            |
| 3              | XX         | XX          | XX            |
| 4              | XX         | XX          | XX            |
| 5              | XX         | XX          | XX            |
| 6              | XX         | XX          | XX            |
| 7              | XX         | XX          | XX            |
| 8              | XX         | XX          | XX            |

## 4 Running MISMIP+ Model 1

### 4.1 Configured settings
To confirm all settings are correctly configured in your environment and you are passing correct arguments to the execution script, Step 0 prints a log of key configuration options. Before proceeding with the model parameterisation and execution, it is recommended to run Step 0, as follows:

```bash
cd ~/ACCESS-ISSM/examples/mismip/
python run-mismip.py --project_code <PROJECT_CODE> --steps 0 --execution_dir <EXECUTION_DIR> --storage <STORAGE_LOC>
```

Without changing any of the default options, this should generate an output similar to the following:

```bash
=============================================================
 ACCESS-ISSM MISMIP+ CONFIGURATION SETTINGS 
=============================================================
 ISSM directory:              /g/data/vk83/prerelease/apps/spack/0.22/release/linux-rocky8-x86_64/gcc-13.2.0/issm-git.7f75a90562fae6834a33bdc6276eaa5dc2cbcf66_0-git.1990-zn2sczssqy2vfrzcdww6ctl5ziinasap
 Python executable:           /g/data/xp65/public/apps/med_conda_scripts/analysis3-25.09.d/bin/python3
---------------------------------------------
 Project code:                <PROJECT_CODE>
 User login:                  <LOGIN>
 Storage locations:           <STOARAGE_LOC>+gdata/vk83
 Model number:                1
 Execution directory:         <EXECUTION_DIR>
 Model name:                  mismip_model_1
 Output directory:            <EXECUTION_DIR>/mismip_model_1
 Steps to run:                [0]
 Load only:                   False
 Walltime:                    2880
 Queue:                       normal
=============================================================
```

### 4.1 Model parameterisation
Steps 1 and 2 are both required for complete model parameterisation. These steps can typically be executed on a _Gadi_ login node and do not require a PBS job submission. To run steps 1 and 2, run:

```bash
cd ~/ACCESS-ISSM/examples/mismip/
python run-mismip.py --project_code <PROJECT_CODE> --steps 1 2 --execution_dir <EXECUTION_DIR> --storage <STORAGE_LOC>
```

This should generate two NetCDF files in `$EXECUTION_DIR/mismip_model_1/`:
* `mismip_model_1_mesh.nc ` conatins the mesh information generated in Step 1.
* `mismip_model_1_parameterise.nc` contains the mesh information generated in Step 1 and the parameterised fields from Step 2.