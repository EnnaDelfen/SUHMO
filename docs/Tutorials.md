
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


# Tutorials: SHMIP test cases

The Subglacial Hydrology Model Intercomparison Project[^1], or SHMIP, is built around six synthetic Suites of experiments (labelled from A to F), each consisting of a set of four to six numerical experiments, designed to show the formation and evolution of the different subglacial drainage elements (sheets and channels) in the context of different input scenarios. Two geometries are considered, a land-terminating ice sheet margin and a synthetic valley-glacier geometry, as shown in Fig a) and b), respectively.

<p align="center">
  <img width="500" src="http://ennadelfen.github.io/SUHMO/IMG/SHMIPgeometries.png">
</p>

This section will walk you trough setting up and running SUHMO's version of the SHMIP suite of test cases A and B. 


## Reminder: setting up your dependencies and environment

Make sure to visit the [Quick start page](https://ennadelfen.github.io/SUHMO/GettingStarted) ! You will find all information relative to setting up your environment, including downloading all dependencies. 


## Suite A

Suite A is a distributed-flow test case, with a steady and spatially uniform water input. The input increases by four orders of magnitude from a low value corresponding to basal melt production (Run A1, m ≃ 2.5 mm a$^{-1}$) to a high water input.

The topography considered is the synthetic representation of a land-terminating ice sheet margin (see a) in the previous image). The ice-sheet domain measures 100 km in the $x$ direction and 20 km in the $y$ direction, the bed is flat and a parabolic ice surface $z_s$ is prescribed by:

$$ z_s(x,y) =  6 (\sqrt{x + 5000} - \sqrt{5000} ) + 1 $$

All you have to do to run the 6 numerical experiments of Suite A is go into the exec folder `A_SHMIP` and build the executable:

```
cd $SUHMO_HOME/exec/A_SHMIP
make all
```

Then you can go into any of the subfolders `AX`, where X=1...6 and launch the computation(s):

```
cd AX
../Suhmo2d.**.ex  input.hydro
```

All simulations are single level, and use a $dt$ of 1 hour, and a $\Delta_x$ of 312.5 m. Results for each numerical experiment, following the format of the SHMIP Supplementary material[^2], is provided in the following sections. 

### A1

<p align="center">
  <img width="800" src="http://ennadelfen.github.io/SUHMO/IMG/SHMIPA1.png">
</p>

### A2

<p align="center">
  <img width="800" src="http://ennadelfen.github.io/SUHMO/IMG/SHMIPA2.png">
</p>

### A3

<p align="center">
  <img width="800" src="http://ennadelfen.github.io/SUHMO/IMG/SHMIPA3.png">
</p>

### A4

<p align="center">
  <img width="800" src="http://ennadelfen.github.io/SUHMO/IMG/SHMIPA4.png">
</p>

### A5

<p align="center">
  <img width="800" src="http://ennadelfen.github.io/SUHMO/IMG/SHMIPA5.png">
</p>

### A6

<p align="center">
  <img width="800" src="http://ennadelfen.github.io/SUHMO/IMG/SHMIPA6.png">
</p>


## Suite B

Comming soon ...

### B1
### B2
### B3
### B4
### B5
### B6



[^1]:de Fleurian, B., Werder, M. A., Beyer, S., Brinkerhoff, D. J., Delaney, I., Dow, C. F., Downs, J., Gagliardini, O., Hoffman, M. J., Hooke, R. L., et al.: SHMIP The subglacial hydrology model intercomparison Project, Journal of Glaciology, 64, 897–916, [doi](https://doi.org/10.1017/jog.2018.78), 2018.
[^2]:de Fleurian, B., Werder, M. A., Beyer, S., Brinkerhoff, D. J., Delaney, I., Dow, C. F., Downs, J., Gagliardini, O., Hoffman, M. J., Hooke, R. L., et al.:Supplementary Material, [doi](https://doi.org/10.1017/jog.2018.78)
