# SUHMO


## Overview

SUHMO is an adaptive-mesh code for simulating the 2D evolution of subglacial hydrology. SUHMO is parallelized with MPI for CPUs. 
Support for complex geometries so far is minimal.

SUHMO relies on the Software for Adaptive Solutions of Partial Differential Equations CHOMBO 
([home page](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)).

## Documentation

A partial documentation with tutorials for SUHMO is available [here](https://ennadelfen.github.io/SUHMO/).

### Getting started

A first simple 2D channelizing problem is available under the SUHMO QuickStart section:

https://ennadelfen.github.io/SUHMO/GettingStarted.html

### Core Algorithm

The SUHMO governing equations and core algorithm are described in:

https://ennadelfen.github.io/SUHMO/Model.html

### SHMIP test cases

A set of self-contained tutorials describing more complex problems is also provided:

https://ennadelfen.github.io/SUHMO/Tutorials.html

## Contributing

New contributions to SUHMO are welcome !

The SUHMO contributions workflow follows these steps:
1. Fork the main repository
2. Create an `AmazingNewFeature` branch implementing your changes 
3. Open a Pull Request from `AmazingNewFeature` on your fork to branch `development` of the main SUHMO repository

Follow [GitHub directions](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) 
to fork SUHMO main repo on your GitHub account, then clone the SUHMO dependencies 
([CHOMBO](https://github.com/EnnaDelfen/Chombo_3.2/tree/feature_SUHMO)) along with your own *SUHMO* fork on your local machine.

Then step into the SUHMO folder and add the main SUHMO repository as the `upstream` remote in order to keep track of the main repo :

       git add remote upstream https://github.com/EnnaDelfen/SUHMO

At any point, you can update the `developement` branch of your local repository with changes implemented in the main repo by pulling from `upstream` : 

        git checkout development
        git pull upstream development

You are now free to modify your own fork! To add a new feature, the procedure is:

1. Create a branch for the new feature from the `development` branch (locally):

        git checkout development 
        git checkout -b AmazingNewFeature

2. and commit your changes to your local repo: 

        git commit -m "Developed AmazingNewFeature"

3. Alongside your development, regularly merge changes from the main repo `development` branch into your `AmazingNewFeature` branch,
fix any conficts, and push your changes to your GitHub fork:
   
        git pull upstream development     [Fix arising conflicts]
        git push -u origin AmazingNewFeature 

4. When you are ready to propose your new feature/improvement/bug fix to the main SUHMO repo, reiterate Step 3 and submit a Pull Request through the GitHub page from your fork onto the `development` branch of the main repo:

 - Click on the ``compare & pull request`` button to start your PR.
 - Provide a title and a short description for your PR:
   * what feature/fix do you propose
   * how did you test it
   * any other information deemed useful : does it modify the default SUHMO behavior ? Does it break retro compatibility for certain tests ? ...
 - Press ``Create pull request``.

Please DO NOT write large Pull Requests, as they are very difficult and time-consuming to review.
As much as possible, split them into small targeted PRs.
For example, if you find typos in the documentation open a pull request that only fixes typos.
If you want to fix a bug, make a small pull request that only fixes this bug.

## Acknowledgment

