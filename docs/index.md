# Welcome to SUHMO's documentation!

- [Getting Started](https://ennadelfen.github.io/SUHMO/GettingStarted) (Quick start section)
- [Core Algorithm](https://ennadelfen.github.io/SUHMO/Model)
- [Tutorials](https://ennadelfen.github.io/SUHMO/Tutorials)

## Overview
SUHMO is an adaptive-mesh code for simulating the 2D evolution of subglacial hydrology. SUHMO is parallelized with MPI for CPUs. Support for complex geometries so far is minimal.
More information about the core algorithm and the governing equations solved can be found on the [dedicated page](https://ennadelfen.github.io/SUHMO/Model).


## Dependencies
SUHMO relies on the Software for Adaptive Solutions of Partial Differential Equations CHOMBO ([see home page](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)).
Before starting, make sure that all the dependencies are properly installed and environment variables properly set ! This is explained in the [Quick start](https://ennadelfen.github.io/SUHMO/GettingStarted) section, but the main steps are also recalled below:

### Get Chombo v3.2
Start by downloading the `mergeHCSUHMO` branch of the forked version of Chombo 3.2:

```
git clone https://github.com/EnnaDelfen/Chombo_3.2.git
git checkout -b mergeHCSUHMO remotes/origin/mergeHCSUHMO
```

Set up the CHOMBO environment variable:
```
export CHOMBO_HOME=path_to_your_CHOMBO3.2_folder/lib/
```


### Get SUHMO v1.0
Then, get the `master` branch of SUHMO:

```
git clone https://github.com/EnnaDelfen/SUHMO.git
```

Set up the SUHMO environment variable:
```
export SUHMO_HOME=path_to_your_SUHMO1.0_folder/
```


### Additional dependencies
To be successful at building and running SUHMO, you will also need MPI and HDF5 installed somewhere on your system.


## Get Started !

Now that your environment is ready, you can start with a very [basic tutorial](https://ennadelfen.github.io/SUHMO/GettingStarted)
