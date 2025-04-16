# ACCESS-ISSM

## About the model
The **ACCESS Ice Sheet System Model (ACCESS-ISSM)** is an implementation of NASA’s [Ice Sheet System Model (ISSM)](https://issm.jpl.nasa.gov/). It is designed for high-resolution, high-performance simulations of ice sheet dynamics, including glacier and ice shelf flow. By leveraging ISSM’s powerful capabilities, ACCESS-ISSM facilitates advanced research into ice sheet behavior, contributing to improved understanding of glaciological processes and their impact on the broader climate system.

For more information see the [ACCESS-Hive Docs model description](#) and how to run the model.

## About this repository
This is the **Model Deployment Repository** for the **ACCESS-ISSM** model.

## Releases
Release information is available on the [ACCESS Hive Forum release topic](#).

## Support
ACCESS-NRI supports ACCESS-ISSM for the Australian Research Community.

Any questions about ACCESS-NRI releases of ACCESS-ISSM should be posted on the [ACCESS-Hive Forum](#). See the [ACCESS Help and Support topic](#) for details on how to do this.

## Build
ACCESS-NRI uses **Spack**, a build-from-source package manager designed for use with high-performance computing. This repository contains a `spack.yaml` environment file that defines all the essential components of the model, including exact versions.

Spack automatically builds all the components and their dependencies, producing model component executables. It already contains support for compiling thousands of common software packages. Spack packages for the components are defined in the [spack packages repository](#).

ACCESS-ISSM is built and deployed automatically to **gadi** on **NCI** (see below). However, it is possible to use Spack to compile the model using the `spack.yaml` environment file in this repository. To do so, follow the instructions for [configuring Spack on gadi](#).

Then clone this repository and run the following commands on **gadi**:

```bash
spack env create access-issm spack.yaml
spack env activate access-issm
spack install
```

This creates a Spack environment called `access-issm` and builds all the components. Their locations can be found by running:

```bash
spack find --paths
```

## Deployment
ACCESS-ISSM is deployed automatically when a new version of the `spack.yaml` file is committed to the `main` branch or a dedicated `backport/VERSION` branch. All the ACCESS-ISSM components are built using Spack on **gadi** and installed under the `vk83` project in `/g/data/vk83`. It is necessary to be a member of the `vk83` project to use ACCESS-NRI deployments of ACCESS-ISSM.

The deployment process also creates a GitHub release with the same tag. All releases are available under this repository’s **Releases** page. Each release has a changelog and metadata with detailed information about the build and deployment, including:

- Paths on **gadi** to all executables built in the deployment process (`spack.location`)
- A `spack.lock` file, which is a complete build provenance document listing all the components that were built and their dependencies, versions, compiler version, build flags, and build architecture. It is also installable via Spack, similarly to the `spack.yaml`.
- The environment `spack.yaml` file used for deployment

Additionally, the deployment creates environment modulefiles, the standard method for deploying software on **gadi**. To view available ACCESS-ISSM versions:

```bash
module use /g/data/vk83/modules
module avail access-issm
```

For users of ACCESS-ISSM model configurations released by ACCESS-NRI, the exact location of the model executables is not required. Model configurations will be updated with new model components when necessary.

