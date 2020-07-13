#!/usr/bin/env bash

#>  Copyright (c) 2020
#>  Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#>  Schloss-Wolfsbrunnenweg 35
#>  69118 Heidelberg, Germany
#
#>  Please send your contact address to get information on updates and
#>  new features to "mcmsoft@h-its.org". Questions will be
#>  answered as soon as possible.
#>  References:
#>  A rapid identification of hit molecules for target proteins via physico-chemical descriptors.
#>  (2013) Phys. Chem. Chem. Phys., 15, 9107-9116.
#>  Authors: Goutam Mukherjee and B. Jayaram
#>  Version 1.0 (April 2013)

#>  RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features
#>  Authors: Stefan Holderbach, Lukas Adam, B. Jayaram, Rebecca C. Wade, Goutam Mukherjee
#>  Version 1.0 (June 2020)
#>  ChemRxiv preprint (https://doi.org/10.26434/chemrxiv.12636704.v1), 2020

# Run this script once to compile and install the dependencies for RASPD+ from this directory
if [ ! -e raspdml.yml ]
then
	echo "Please cd to the root of the RASPD+ repository for this installation"
	exit 1
fi

if [ -z "$raspd_root" ]
then
    source config/init.sh
fi

if [ -z "$raspd_root" ] 
then
	echo "Please set raspd_root to the root of your RASPD+ installation"
	echo "   export raspd_root=<where RASPD+ is installed>"
	exit 1
fi

if [ -z "$conda_root" ] 
then
	echo "Please set conda_root to the root of your conda installation"
	echo "   export conda_root=<where conda is installed (..../miniconda)>"
	exit 1
fi

eval "$($conda_root/condabin/conda shell.bash hook)"

conda env create -f $raspd_root/raspdml.yml

conda env create -f $raspd_root/trappenv.yml

mkdir -p $raspd_root/bin

cd $raspd_root/src

make

cd $raspd_root

wget 'https://zenodo.org/record/3937426/files/weights.tar.gz'
tar -xzf weights.tar.gz

