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
#>  Manuscript in preparation, 2020

# Code to manage data splitting and preprocessing for RASPD+
# Author(s): Stefan Holderbach, Lukas Adam
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import os
import numpy as np
from sklearn.preprocessing import RobustScaler
from sklearn.model_selection import train_test_split, KFold, ShuffleSplit


class DataPrep:

    def __init__(self, input, test_size=0.1, n_splits=4, shuffle=True, single_mode=True, n_shuffle=1, random_state=42):
        """
        Helper class to perform splitting of datasets into partitions
        :param input: Tuple of (x, y, index) with sample in first dimension
        :param test_size: Ratio of data used for hold out test set
        :param n_splits: number of inner k-fold splits for train-val CV
        :param shuffle: Has to be true if you want multiple test-train splits
        :param single_mode: Only perform one train test split (Affects output dimensions!)
        :param n_shuffle: Number of shuffled splits for test set vs train-val data (only used if single_mode=False)
        :param random_state: Seed for shuffling and splitting
        """

        if len(input) != 3:
            raise ValueError('no data was forwarded to split. Please input x, y and indices')
        else:
            self.x, self.y, self.index = input

        self.test_size = test_size
        self.n_splits = n_splits
        self.shuffle = shuffle
        self.random_state = random_state

        self.single_mode = single_mode
        self.n_shuffle = n_shuffle if not single_mode else 1
        self.scalers = [RobustScaler() for _ in range(self.n_shuffle)]
        self.x, self.x_test, self.y, self.y_test, self.index, self.index_test = self.split_test_rest()
        self.scale_data(self.x, mode='fit')
        self.x = self.scale_data(self.x, mode="transform")
        self.x_test = self.scale_data(self.x_test, mode='transform')

        self.train_val_test = {
            'x_train': [],
            'x_val': [],
            'y_train': [],
            'y_val': [],
            'index_train': [],
            'index_val': [],
        }
        self.make_train_val_test()

    def split_test_rest(self):
        """
        Splits the dataset into train and test
        Takes self.x, self.y, self.index
        :return: x, x_test, y, y_test, index, index_test
        """
        # Ensures that equally sized folds result
        train_size = round((1 - self.test_size) * self.x.shape[0])
        train_size -= train_size % self.n_splits
        test_size = self.x.shape[0] - train_size

        if self.single_mode:
            split_test_rest = train_test_split(self.x, self.y, self.index, test_size=test_size,
                                               random_state=self.random_state)
        else:
            jf = ShuffleSplit(n_splits=self.n_shuffle, test_size=test_size, random_state=self.random_state)
            indices = [*jf.split(self.x, self.y)]
            split_test_rest = [[] for _ in range(6)]
            for i, dat in enumerate((self.x, self.y, self.index)):
                for train, test in indices:
                    split_test_rest[i * 2].append(dat[train])
                    split_test_rest[i * 2 + 1].append(dat[test])
        return split_test_rest

    def scale_data(self, x, index=0, mode='fit'):
        """
        Helper function for scaling data, handling scalers.
        Works on a single data array or a list of data arrays with the same length as self.n_shuffle
        :param x: Data shape [(n_shuffle,) n_samples, n_features]
        :param index: Which of the scalers to use (ignored if given array/list with n_shuffle in the first dim)
        :param mode: {"fit", "transform", "reverse"}
        :return: None or (list of) rescaled data
        """
        if self.single_mode:
            assert index == 0, "Index can only be 0 for single mode"

        x_rank = np.array(x).ndim
        if x_rank == 3:
            assert len(x) == self.n_shuffle, "Data with same number of replicates as n_shuffle expected"
            if mode == 'fit':
                for i in range(self.n_shuffle):
                    self.scalers[i].fit(x[i])
            elif mode == 'transform':
                return [self.scalers[i].transform(x[i]) for i in range(self.n_shuffle)]
            elif mode == "reverse":
                return [self.scalers[i].inverse_transform(x[i]) for i in range(self.n_shuffle)]
            else:
                raise ValueError('method for scaling not defined')
        elif x_rank == 2:
            if mode == 'fit':
                self.scalers[index].fit(x)
            elif mode == 'transform':
                return self.scalers[index].transform(x)
            elif mode == "reverse":
                return self.scalers[index].inverse_transform(x)
            else:
                raise ValueError('method for scaling not defined')
        else:
            raise ValueError("x either has to be of rank 1, 2 with a specified index or of rank 3")

    def make_train_val_test(self):
        kf = KFold(n_splits=self.n_splits, shuffle=self.shuffle)
        data = {
            'x': self.x,
            'y': self.y,
            'index': self.index
        }
        if self.single_mode:
            train_indxs, val_indxs = zip(*kf.split(self.x))
            for index, elem in data.items():
                to_stack = [elem[idxs] for idxs in train_indxs]
                self.train_val_test[index + '_train'] = np.stack(to_stack, axis=0)
                to_stack = [elem[idxs] for idxs in val_indxs]
                self.train_val_test[index + '_val'] = np.stack(to_stack, axis=0)
        else:
            for key in self.train_val_test:
                self.train_val_test[key] = [[] for _ in range(self.n_shuffle)]
            for i in range(self.n_shuffle):
                train_indxs, val_indxs = zip(*kf.split(self.x[i]))
                for index, elem in data.items():
                    to_stack = [elem[i][idxs] for idxs in train_indxs]
                    self.train_val_test[index + '_train'][i] = np.stack(to_stack, axis=0)
                    to_stack = [elem[i][idxs] for idxs in val_indxs]
                    self.train_val_test[index + '_val'][i] = np.stack(to_stack, axis=0)
            for key, item in self.train_val_test.items():
                self.train_val_test[key] = np.stack(item)
        self.train_val_test['x_test'] = np.stack(self.x_test)
        self.train_val_test['y_test'] = np.stack(self.y_test)
        self.train_val_test['index_test'] = np.stack(self.index_test)

    def save_datasets(self, save_to_dir):
        """
        Stores self.train_val_test into npy files named by the keys
        :param save_to_dir:
        :return:
        """
        for name_set, set in self.train_val_test.items():
            np.save(os.path.join(save_to_dir, name_set), set)

    def retrieve_datasets(self):
        """
        :return: Dict: train_val_test of npy arrays [(n_shuffle,) (n_splits,) sample, dims],
        (n_shuffle list of) scaler object
        """
        return self.train_val_test, self.scalers
