# MicrobiomeAgeModel2024

_ALPHA VERSION_

Analysis code and data acquisition instructions for the manuscript entitled
"Early life microbial succession in the gut follows common patterns in humans across the globe"
by Bottino, G., Bonham, K. S., Klepac-Ceraj, V. et al.

Work supported by the Wellcome Leap 1kD program.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12822333.svg)](https://doi.org/10.5281/zenodo.12822333)
[![Dryad Badge](https://img.shields.io/badge/Data%20Dryad%3A-10.5061%2Fdryad.dbrv15f9z-blue?color=blue)](https://doi.org/10.5061/dryad.dbrv15f9z)

## System requirements

### Julia

This repository is a julia language (https://julialang.org/) project. In order to run the code, you will need to install the julia language in your system. The code in this repository should work with julia v1.9.3 or higher. There are several ways to install julia, we recommend using `juliaup`, or downloading and installing from the julia website.

Julia is available for Windows, Mac and several Linux distributions out-of-the-box. We tested this code on Microsoft Windows 10 Pro and Linux Mint 20.1 "Ulyssa" LTS.

## Installation Guide

### Project environment and dependencies

Once julia is installed, you can download this repository using `git`, or from the [latest Zenodo relase](https://doi.org/10.5281/zenodo.12822333).
Then, run julia, and [activate this project](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).
If julia is installed in your system `PATH`, you can activate it directly from the command line using `julia --project=@.`. 
For example, if you downloaded and unpacked the code in your `Documents` directory,
you can run:

```
$ cd Documents/MicrobiomeAgeModel2024

$ julia --project=@.
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.3 (2023-08-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

Alternatively, you can change directories from the julia REPL,
and then activate the project from the REPL:

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.3 (2023-08-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> cd(expanduser("~/Documents/MicrobiomeAgeModel2024))

julia> # press ] to enter the Pkg REPL

(@v1.8) pkg> activate .
Activating project at `~/Documents/MicrobiomeAgeModel2024`

(MicrobiomeAgeModel2024) pkg>
```

Either way, once you have the project activated, run `instantiate` from the `Pkg` REPL

```
(MicrobiomeAgeModel2024) pkg> instantiate
  Installed SIMDDualNumbers ────────────────── v0.1.1                                                     
  Installed StatsFuns ──────────────────────── v1.1.1
  Installed CategoricalDistributions ───────── v0.1.9
  Installed JpegTurbo_jll ──────────────────── v2.1.2+0
  #... etc
```

This step only has to be executed once and typically takes less than 5 minutes on a modern laptop or workstation.

### Leap.jl package

One of the project dependencies is the `Leap.jl` package, developed by our lab in the context of the Leap grant to centralize APIs to deal with Machine Learning on biological data. This package is not present in the julia registry, so it has to be added manually from the github repository. To do so, activate the project and enter package mode (press `]`), then type

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.3 (2023-08-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> # press ] to enter the Pkg REPL

(MicrobiomeAgeModel2024) pkg> add https://github.com/Klepac-Ceraj-Lab/Leap
```

## Running the code

We provided a series of markdown notebooks that allow any user to reproduce the code that generated the results, figures and tables in the original manuscript. The main notebooks can be found on the `manuscript` folder of the Repo. These notebooks are divided by numbered figure in the main text (1,2,3 and 4) - with a few additional auxiliary notebooks on the `notebooks` folder (_to be released_). Each notebook has comments to walk the user through the necessary inputs and expected outputs. Outputs will by default be written to a `results` folder in the package root directory.

The execution time for all four notebooks is dominated by the model training with hyperparameter grid search, and I/O operations. Expect the full CV training time to take 1-2 hours on a modern lapop or workstation, and the same time for the combination of LOOCV experiments.