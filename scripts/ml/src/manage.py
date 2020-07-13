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

# Utilities to manage unique folders and a logger instance
# Author(s): Lukas Adam, Stefan Holderbach
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import os
import logging
import datetime


class ManagementSystem:
    # TODO: Declass this convoluted function
    def __init__(self, datadir=None, make_unique=True):
        if datadir is None:
            raise ValueError('no data directory was forwarded. Please set path')
        else:
            self.datadir = datadir

        print('\nStart logger')
        self.logger = logging.getLogger(__name__)

        print('\nStart generating folder')
        self.unique_folder = self.generate_folder() if make_unique else self.datadir
        os.makedirs(self.unique_folder, exist_ok=True)

        print('\nSetting basic configuration for logger to: ',
        os.path.join(self.unique_folder, 'logfile.log'), 'level: ', logging.INFO)
        logging.basicConfig(filename=os.path.join(self.unique_folder, 'logfile.log'), level=logging.INFO)

        if make_unique:
            self.logger.info('Generated unique folder: ' + self.unique_folder)

    def generate_folder(self):
        suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
        basename = "raspdata"
        foldername = "_".join([basename, suffix])
        return os.path.join(self.datadir, foldername)

    def setting(self):
        print('\nget current setting')
        return self.unique_folder, self.logger
