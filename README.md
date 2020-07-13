# RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features

See the preprint on ChemRxiv: https://doi.org/10.26434/chemrxiv.12636704.v1

This repository contains the code used to create the RASPD+ descriptors and machine learning models and shell scripts for approximate protein ligand binding energy useful for prefiltering in virtual screening experiments.

Additional data and the associated model weights can be found on zenodo (https://doi.org/10.5281/zenodo.3937426)

## Documentation

A user documentation documentation can be found at `doc/ReadMe-RASPDPlus.pdf`. For additional technical information check the associated manuscript as well as ReadMe files in the subdirectories.

## Installation

### Requirements

* unixoid system (tested on Ubuntu 18.04 LTS)
* `bash` shell to run the CLI scripts
* `conda` for management of python dependencies (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
* `gcc`, `g++`, `make`
* `TRAPP` >= 4.0.1 (https://www.h-its.org/downloads/trapp/)

### Instructions

Please follow the instructions provided in the user documentation document

## Folder structure

### `src`

This folder contains all C(++)-programs.
These are required to calculate the physicochemical properties of ligand and ligand binding pocket of a protein

### `bin`

This folder contains all executable files (Created by running `install.sh`)

### `scripts`

This folder contains all scripts which call all the executable files from bin folder to calculate
the physicochemical properties of the ligand and the binding pocket of a protein and run the machine learning models to predict the binding free energy.
To estimate the volume of protein pocket it calls TRAPP (TRAnsient Pockets in Proteins).

The folder `ml` contains the code to train and evaluate the machine learning models as well as jupyter notebooks which were used to create the figures in the paper.

### `weights`

In this folder the model weights are placed. A download of the primary model weights from the zenodo store (https://doi.org/10.5281/zenodo.3937426) is automatically performed upon running `install.sh`

### `data`

This folder contains the

(a) physicochemical parameters (`PARAMETER_MW`, `maxdistance` and `PDBID`) of million molecules which were retrieved from ZINC database;

(b) AMBER input file (`s1.cmd`) to add hydrogen atoms to protein

(c) SMILES code for the million molecules (`target.txt`)

(d) a file (`select_parameter.txt`) to configure the acceptable range of physicochemical parameters while screening million or customized data set.

### `doc`

In this folder the RASPD+ user manual can be found.

### `customized_data`

This folder is empty. However, the `lig_parameters_gen.sh` script generates physicochemical properties of ligands which are autosaved in this folder with a fixed file name

(a) `list_customized` (contains molecule ids)

(b) `PARAMETER_MW_Customized` (contains physicochemical parameters of molecules)

(c) `maxdistance_customized`  (contains DMax values of molecules)

(d) `target_customized.smi` (contains SMILES strings of molecules)

(e) `s1.cmd` (AMBER input)

### `example`

This folder contains a few example files, and results of testruns and a ReadMe file how to reproduce them with RASPD+.

### `config`

In this folder a script file, `init.sh`, resides. It is required to source the root of the cloned git repository
and path to your conda installation eg, `/home/<user name>/miniconda3`
Alternatively the configuration provided by this script can be moved to your bash profile (e.g. `~/.bashrc`)

## Copyright

Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org) \
Schloss-Wolfsbrunnenweg 35\
69118 Heidelberg, Germany

Supercomputing Facility for Bioinformatics and Computational Biology, IIT Delhi (http://www.scfbio-iitd.res.in/) \
Indian Institute of Technology\
Hauz Khas, New Delhi - 110016, India

Authors of RASPD+ (version 1.0 (June 2020)): Stefan Holderbach, Lukas Adam, B. Jayaram, Rebecca C. Wade, Goutam Mukherjee

Authors of the original RASPD (version 1.0 (April 2013)): Goutam Mukherjee and B. Jayaram

## License

This software is licensed under the European Union Public License version 1.2. For details see `LICENSE`

## Contact

Please send an email to get information on updates and new features to "mcmsoft@h-its.org".
Questions will be answered as soon as possible.
