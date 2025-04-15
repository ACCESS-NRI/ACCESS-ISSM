# ISSM Model Deployment Repository

This repository contains the **model deployment** configuration for **ISSM** using [Spack](https://spack.readthedocs.io/en/latest/). Below are instructions and considerations for setting up the necessary files and configurations. Please follow these steps carefully to ensure consistency and reproducibility of the deployment environment.

---

## 1. GitHub Actions Configuration (`.github/workflows.yaml` or `config/workflows.yaml`)

If you use GitHub Actions workflows to validate or build your deployment, note the following:

1. **Root SBD Name**  
   If the name of the _root Spack Build Descriptor (SBD)_ for ISSM differs from the repository or model name (e.g., `access-esm1p5` is used for `ACCESS-ESM1.5`), you need to:
   ```yaml
   jobs:
     pr-ci:
       with:
         # root-sbd: <ROOT_SBD_NAME>
     pr-comment:
       with:
         # root-sbd: <ROOT_SBD_NAME>
     pr-closed:
       with:
         # root-sbd: <ROOT_SBD_NAME>
Uncomment the root-sbd line(s).
Set each root-sbd to the correct name of the root SBD you want Spack to install.
Workflow Steps
Ensure that your workflow checks out the relevant files and directories (e.g., config/, spack.yaml) and that it uses the correct environment variables or inputs when calling spack.
2. Versioning Configuration (config/versions.json)

Your config/versions.json file should specify versions for both:
```
.spack
Provide a semver or CalVer version, typically matching a branch in the ACCESS-NRI/spack repository. For example:
{
  ".spack": "1.2.3",
  ...
}
```
This will be used to clone the releases/v1.2.3 branch of ACCESS-NRI/spack, if following typical usage patterns.
.spack-packages
Provide a CalVer-compliant version tag that matches one of the tags in ACCESS-NRI/spack-packages. For example:
```
{
  ".spack": "1.2.3",
  ".spack-packages": "2024.01.01",
  ...
}
```
3. Main Spack Configuration (spack.yaml)

Below are the key sections and how to configure them for ISSM. Adjust these instructions based on whether you need single or multiple compiler/target variants.

3.1 If There Are Multiple Deployment Targets
If you need to build ISSM with different compilers, targets, or other variants across multiple deployments, you can use definitions in spack.yaml:
```
spack:
  definitions:
    - ROOT_PACKAGE: issm
    - compiler_target:
      - gcc@11.2.0
      - nvhpc@21.9

  specs:
    - $ROOT_PACKAGE  # This references spack.definitions above.
                     # At build time, $ROOT_SPEC might be formed using 
                     # the appropriate compiler/target combos.
```
  # More configuration below...
Typically, you will construct a ROOT_SPEC from $ROOT_PACKAGE + variants + compiler, etc. Each environment matrix entry can expand to its own spec. Adjust this pattern as necessary.

3.2 If There Is Only a Single Compiler/Target
If you do not require multiple compilers or target variations, simplify spack.yaml:
```
spack:
  specs:
    - issm@git.2025.04.0

```
issm should be the root SBD name (the package name in spack-packages).
@git.2025.04.0 follows a CalVer-like scheme to indicate a specific tag or version for this entire deployment.
3.3 Packages and Variants
In the spack.packages section, you can specify dependencies, versions, and variants. For example:
```
spack:
  packages:
    netcdf:
      require:
        - 4.8.1
      variants:
        - +mpi
        - +parallel-netcdf
    hdf5:
      require:
        - 1.12.1
      variants:
        - +mpi
    ...
```
The first element of each require list must be only the version number. Any variants or additional constraints go under variants (or separate attributes).
3.4 All-Package Configuration
You can also set global configuration under spack.packages.all, such as:
```
spack:
  packages:
    all:
      compiler: [gcc@11.2.0]
      target: [x86_64]
```
3.5 Module Configuration
To make certain packages available in module form (`module load <package>`), edit the spack.modules.default.tcl section. For example:
```
spack:
  modules:
    default:
      tcl:
        include:
          - issm  # Packages you want to be directly loadable
          - netcdf
          - hdf5

        projections:
          "issm": "{name}/{version}"
          "netcdf": "{name}/{version}"
          "hdf5": "{name}/{version}"
```
For each included module, you must set the name of the module to match spack.packages.*.require[0] (the version) – but typically in the simplest form (e.g., 4.8.1 for NetCDF).

4. Repository Notes

Ensure this repository has the correct version tags (e.g., git.2025.04.0) that match the @git.YEAR.MONTH.MINOR references in spack.yaml.
If you rename the root SBD or require a different package name for ISSM, update:
config/workflows.yaml (or your GitHub Actions) for root-sbd.
spack.yaml for the specs section and any references in spack.packages.
5. Contributing

Please submit issues and pull requests if you encounter problems or have improvements. For major changes, open an issue first to discuss proposed modifications.

6. License

This repository is provided under the MIT License (or whichever license applies). Refer to the LICENSE file for details.

Happy deploying ISSM with Spack!
Feel free to reach out if you have any questions or suggestions.