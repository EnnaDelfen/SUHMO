
<head>

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.0/dist/katex.min.css" integrity="sha384-Xi8rHCmBmhbuyyhbI88391ZKP2dmfnOl4rT9ZfRI7mLTdk1wblIUnrIq35nqwEvC" crossorigin="anonymous">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.0/dist/katex.min.js" integrity="sha384-X/XCfMm41VSsqRNQgDerQczD69XqmjOOOwYQvr/uuC+j4OPoNhVgjdGFwhvN02Ja" crossorigin="anonymous"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.0/dist/contrib/auto-render.min.js" integrity="sha384-+XBljXPPiv+OzfbB3cVmLHf4hdUFHlWNZN5spNQ7rmHTXpd7WvJum6fIACpNNfIR" crossorigin="anonymous"></script>
<script>
    document.addEventListener("DOMContentLoaded", function() {
        renderMathInElement(document.body, {
          // customised options
          // • auto-render specific keys, e.g.:
          delimiters: [
              {left: '$$', right: '$$', display: true},
              {left: '$', right: '$', display: false},
              {left: '\\(', right: '\\)', display: false},
              {left: '\\[', right: '\\]', display: true}
          ],
          // • rendering keys, e.g.:
          throwOnError : false
        });
    });
</script>
  
</head>


Go back to the [main documentation page](https://ennadelfen.github.io/SUHMO/index)




# SUHMO's quick start section

This section will walk you trough downloading all dependencies and compiling the sources to run a very basic 2D channelizing test case. 

## Dependencies
SUHMO relies on the Software for Adaptive Solutions of Partial Differential Equations CHOMBO ([see home page](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)). Make sure that you create two separate folders, one for Chombo and another one for SUHMO, as it will make following the instructions easier.

### Getting Chombo v3.2
Start by downloading the `feature_SUHMO` branch of the forked version of Chombo 3.2:

```
git clone https://github.com/EnnaDelfen/Chombo_3.2.git
git checkout -b feature_SUHMO remotes/origin/feature_SUHMO
```

Set up the CHOMBO environment variable:
```
export CHOMBO_HOME=path_to_your_CHOMBO3.2_folder/lib/
```

### Getting SUHMO v1.0
Then, get the `master` branch of SUHMO:

```
git clone https://github.com/EnnaDelfen/SUHMO.git
```

Set up the SUHMO environment variable:
```
export SUHMO_HOME=path_to_your_SUHMO1.0_folder/
```

### Additional dependencies
On Ubuntu, use the `apt-get` command to get basic compiler and support libraries (MPI/HDF5):

```
sudo apt-get install -y --no-install-recommends \
    build-essential csh \
    g++-11 gfortran-11 \
    libopenmpi-dev \
    openmpi-bin \
    libhdf5-openmpi-dev
```

Note that some of these might already be installed on your system, so adapt this command line accordingly.


### Adapting the Make.defs.local
Go into the Chombo make folder:

```
cd $CHOMBO_HOME/mk
```

and modify the local `local/Make.defs.template` to provide the requested paths, according to your installation. Generate a symbolik link pointing to it as follows:

```
rm -rf Make.defs.local
ln -s local/Make.def.template Make.defs.local
```

You should be all set to compile and run your first example ! 


## Setting up the channelizing test case

### Test case description

<p align="center">
  <img width="500" src="http://ennadelfen.github.io/SUHMO/IMG/ConvergenceChannelTC.png">
</p>

The domain is a rectangle of 64 m in the $x$ direction by 16 m in the $y$ direction. The bed is sloped in the $x$ direction (with a +2$\%$ slope) and is topped with a slab of ice of constant thickness everywhere (500 m). A moulin delivering a total input of 30 m$^3$ s$^{-1}$ of water at steady state is located 16 m from the outlet of the domain in the $x$ direction. The geometry is shown above. The moulin influx is ramped up over a period of a little under a month. 
Periodicity is assumed in the $y$ direction, a homogeneous Neumann boundary condition (outlet) is prescribed at x$=$0 m while a homogeneous Dirichlet boundary condition is prescribed at x$=$64 m. Note that the moulin source term follows a Gaussian distribution.


### Building the executable
Start by going into the following SUHMO test case folder and building the executable:

```
cd $SUHMO_HOME/exec/0_convergence_channelized
make all
```

### Running the first example
Go into the `1lev` test case and launch your first SUHMO simulation:

```
cd 1lev
../Suhmo2d.**.ex  input.hydro
```

If everything goes according to plan you should not have to wait more than a minute for the simulation to complete. The following files should have been generated:

- A `pout.0` output file filled with run information
- A check file `chk003000.2d.hdf5` to restart your simulation if needed
- 2 plot files `plot00XXXX.2d.hdf5` corresponding to the initial and final states of your simulation

### Analysing the results
There are several tools that enable you to visualize simulations results. [ParaView](https://www.paraview.org/) is an option, as well as [VisIt](https://visit-dav.github.io/visit-website/).


The fields of gap height and head resulting from this (coarse) simulation are shown in the top row of the following image. Note that the discretization is 32x8, leading to a $\Delta_x$ of 2 m. This test case can be used to run a convergence study, and results from a much finer mesh are also shown, in the bottom row. In that case, the discretization is 512x128 so that $\Delta_x$ is 0.125 m (input files for that run are provided in the `5lev` test case folder). 

<p align="center">
  <img width="800" src="http://ennadelfen.github.io/SUHMO/IMG/Solution1lev_5lev.png">
</p>

The `pout.0` file can also be checked for useful information. The number of Picard iterations, FASMG iterations (for $h$) and MG iterations (for $b$) are reported for each time step. The first 50 timesteps exhibit from 2 to 3 Picard interations and over 30 FASMG iterations while the initial state gets settled. Then the moulin input ramps up and as many as 7 Picard iterations are required for another 200-300 iterations, while the channel develops. Steady state is reached soon after.


### Go a bit further
Congrats on running your first SUHMO test case ! 

You will note many subfolders in the `0_convergence_channelized` exec folder. All subfolders labelled `Xlev`, where X is a number, are single level simulations. Subfolders labelled `Xlev_base` are 2 level simulations, where the base level refinement is as for the `Xlev`. Subfolders labelled `Xlev_base2levs` are 3 level simulations, where the base level refinement is as for the `Xlev`.
