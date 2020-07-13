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

# Dataloader for descriptor tables used in RASPD+
# Author(s): Lukas Adam, Stefan Holderbach
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import os
import glob
import numpy as np
import pandas as pd
from sklearn.utils import shuffle


class DataLoader:
    # TODO: Declass this convoluted function
    # TODO: Extract the fix column names as parameters
    def __init__(self, datadir=None, filenames=None,
                 logger=None, to_shuffle=False, read_in_with="txt", ignore_cols=None):

        if logger is None:
            raise ValueError('no logger specified. Please set logger!')
        else:
            self.logger = logger

        if datadir is None:
            raise ValueError('no data directory was forwarded. Please set path')
        else:
            self.datadir = datadir

        if filenames is None:
            self.filenames = glob.glob(os.path.join(self.datadir, '*' + read_in_with))
        elif type(filenames) is list and len(filenames) >= 1:
            self.filenames = filenames
        else:
            raise ValueError('no filenames were specified')

        self.ignore_cols = ['PA(LYN)']  # TODO: Remove hard coding
        if ignore_cols:
            self.ignore_cols.extend(ignore_cols)

        self.data = {}
        for filename in self.filenames:
            try:
                self.logger.info('try to load ' + filename)
                if to_shuffle:
                    self.data[filename] = shuffle(self.read_file(filename))
                else:
                    self.data[filename] = self.read_file(filename)
                self.logger.info('loading ' + filename + ' was succesful')
            except Exception as e:
                print(e)
                continue

    def read_file(self, filename):
        if filename.endswith('.txt'):
            file = pd.read_csv(filename, sep=";")
        elif filename.endswith('.csv'):
            file = pd.read_csv(filename, sep=",")
        elif filename.endswith('.xlsx'):
            file = pd.read_excel(filename)
        else:
            raise ValueError("Unknown filetype for :" + filename)
        # Attention: Hardcoded hotfix
        # TODO: maybe move this line into more RASPD-ML specific place
        file.rename(columns={"Atom_Efficiency": "Atom Efficiency"}, inplace=True)  # Fix the name of atom efficiency
        # Bring the dataset into proper form
        file = file.set_index("PDBID")  # We expect a column named PDBID for labelling
        file = file.astype(np.float)
        file.dropna(inplace=True)  # We ignore bad columns
        assert file.shape[0], filename + " doesn't contain data!"
        return file.drop(columns=self.ignore_cols, errors="ignore")

    def retrieve_data(self):
        return self.data
