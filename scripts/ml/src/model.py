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

# Wrappper code for deep learning model hyperparameter optimization in nested k-fold CV
# Author(s): Lukas Adam, Stefan Holderbach
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import os
import random
import string
import warnings
import keras
import numpy as np
import talos as ta
from talos.utils.gpu_utils import parallel_gpu_jobs

parallel_gpu_jobs()
from scipy import stats
from keras import Sequential
from keras.layers import Dense, Dropout
from keras.models import model_from_json
from keras.optimizers import SGD, Adam, Adadelta, Adagrad, Adamax, Nadam, RMSprop
from talos.model.normalizers import lr_normalizer

warnings.filterwarnings('ignore')


def randomString(stringLength):
    """Generate a random string with the combination of lowercase and uppercase letters """
    letters = string.ascii_letters
    return ''.join(random.choice(letters) for i in range(stringLength))


class ModelKCross:

    def __init__(self, scan, index, y_true, x_pred, fold):
        self.scan = scan
        self.model = None
        self.index = index
        self.y_true = y_true.squeeze()
        self.y_pred = None
        self.x_pred = x_pred
        self.fold = fold
        self.get_model()
        self.run_model()

    def get_saved_model(self):
        return self.scan.saved_models[self.index]

    def get_saved_weights(self):
        return self.scan.saved_weights[self.index]

    def get_model(self):
        self.model = model_from_json(self.get_saved_model())
        self.model.set_weights(self.get_saved_weights())

    def run_model(self):
        self.y_pred = self.model.predict(self.x_pred).squeeze()

    def score(self):
        score, _ = stats.pearsonr(self.y_true, self.y_pred)
        return score


class BestModelKCross:

    def __init__(self, scans):
        self.scans = scans
        self.index = self._index_best_model()
        self.params = scans[0].data.loc[self.index]

    def _index_best_model(self):
        td = []
        for scan in self.scans:
            td.append(np.array(scan.data.loc[:, 'r2_score_val']).astype(float))
        td = np.stack(td, axis=0)
        tb = np.mean(td, axis=0)
        return tb.argmax()

    def get_saved_model(self, scan_index):
        return self.scans[scan_index].saved_models[self.index]

    def get_saved_models(self):
        models = {}
        for scan_index in range(len(self.scans)):
            model = model_from_json(self.get_saved_model(scan_index))
            model.set_weights(self.get_saved_weights(scan_index))
            models['RASPDeep_fold_' + str(scan_index)] = self._compile_model(model)
        return models

    def get_saved_weights(self, scan_index):
        return self.scans[scan_index].saved_weights[self.index]

    def _compile_model(self, model):
        model.compile(loss=self.params['loss_fn'],
                      # here we add a regulizer normalization function from Talos
                      optimizer=make_opt_instance(self.params["opt"], self.params["lr"]))
        return model


def make_opt_instance(opt, lr_scale):
    # Turn the optimizer into a class if necessary
    opt_trans = {
        'sgd': SGD,
        'rmsprop': RMSprop,
        'adagrad': Adagrad,
        'adadelta': Adadelta,
        'adam': Adam,
        'adamax': Adamax,
        'nadam': Nadam,
    }
    _opt = opt_trans[opt.lower()] if isinstance(opt, str) else opt
    return _opt(lr=lr_normalizer(float(lr_scale), _opt))


class DL:

    def __init__(self, x_train, y_train, x_val, y_val, params, wdir):
        self.params = params
        self.x_train = x_train
        self.y_train = y_train
        self.x_val = x_val
        self.y_val = y_val
        self.histories = {}
        self.kfold = x_train.shape[0]
        self.model_index = 0
        self.wdir = wdir
        self.boolean = False
        self.model_name = ''

    def model(self, x_train, y_train, x_val, y_val, params):

        if self.boolean:
            self.model_index = 1
            self.model_name = os.path.join(self.wdir, 'model' + str(self.model_index))
        else:
            self.model_index += 1
            self.model_name = os.path.join(self.wdir, 'model' + str(self.model_index))

        # next we can build the model exactly like we would normally do it
        model = Sequential()
        model.add(Dense(params['node1'],
                        input_shape=(x_train[0].shape[-1],),
                        activation=params['activation1'],
                        kernel_initializer=params['initializer'],
                        bias_initializer='zeros'))
        model.add(Dropout(params['dropout1']))
        model.add(Dense(params['node2'],
                        activation=params['activation2'],
                        kernel_initializer=params['initializer'],
                        bias_initializer='zeros')
                  )
        model.add(Dropout(params['dropout2']))
        # Here you could add more layers and change the parameter dict
        # Compare nn_param_space in train_and_eval.py
        model.add(Dense(1, kernel_initializer=params['initializer'],
                        bias_initializer='zeros'))

        model.compile(loss=params['loss_fn'],
                      optimizer=make_opt_instance(params["opt"], params["lr"]))

        # IF early stopping unreliable check parameters https://keras.io/callbacks/#earlystopping
        earlystopping = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0, patience=30, verbose=0,
                                                      mode='auto')
        callbacks = [earlystopping]
        history = model.fit(x_train, y_train,
                            validation_data=[x_val, y_val],
                            batch_size=params['batch_size'],
                            epochs=params['epochs'],
                            verbose=0,
                            callbacks=callbacks)

        # Get the dictionary containing each metric and the loss for each epoch
        if self.fold in self.histories:
            self.histories[self.fold].append(history.history)
        else:
            self.histories[self.fold] = [history.history]

        self.boolean = False
        # finally we have to make sure that history object and model are returned
        return history, model

    def run_talosScan(self, fold, last_epoch_value=True, experiment_no="1"):
        # and run the experiment
        self.boolean = True
        t = ta.Scan(x=self.x_train[fold],
                    y=self.y_train[fold],
                    model=self.model,
                    params=self.params,
                    x_val=self.x_val[fold],
                    y_val=self.y_val[fold],
                    last_epoch_value=last_epoch_value,
                    experiment_no=experiment_no,
                    search_method='linear')
        return t

    def run(self, max_folds=None):
        self.taloslib = {}
        if max_folds is None:
            max_folds = self.kfold
        elif not 0 < max_folds <= self.kfold:
            max_folds = self.kfold
        for fold in range(max_folds):
            self.fold = str(fold)
            talosScan = self.run_talosScan(fold)
            if fold == 0:
                self.subplot_count = len(list(talosScan.data.index))
            for i in list(talosScan.data.index):
                talosScan.data.loc[i, 'r2_score_val'] = ModelKCross(talosScan, i, self.y_val[fold],
                                                                    self.x_val[fold], fold).score()
            self.taloslib[str(fold)] = {'talosScan': talosScan, 'history': self.histories[self.fold]}

        return self.taloslib

