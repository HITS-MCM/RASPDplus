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

# RASPD+ inference script to run predictions from a saved models folder
# Author(s): Stefan Holderbach
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.model_selection import ParameterGrid
from scipy import stats
from time import perf_counter

DATASET_MODES = ["validation", "test", "training"]


def q2_f3_score(y_test, y_pred, y_train):
    """
    Calculates the $Q^2_{F3} metric for dataset corrected performance evaluation

    :param y_test:
    :param y_pred:
    :param y_train:
    :return:
    """
    return 1 - (metrics.mean_squared_error(y_test, y_pred) / np.var(y_train))


class MLA:
    """
    Machine learning assessment class: performs validation on the different nested datasets for a set of models
    """
    def __init__(self, x_train, y_train, x_val, y_val, x_test, y_test,
                 additional_tests=None, scaler=None, n_perm=0, feature_names=None):
        self.x_train = x_train
        self.y_train = y_train
        self.x_val = x_val
        self.y_val = y_val
        self.x_test = x_test
        self.y_test = y_test
        self.kfold = x_train.shape[0]
        self.score_dict = {}
        self.perm_metrics = {}
        self.add_test = additional_tests
        self.scaler = scaler
        self.n_perm = n_perm
        self.n_features = x_train[0].shape[-1]
        self.feature_names = feature_names

    def run(self, classifier_dict):

        for model_name, models in classifier_dict.items():
            for fold in range(self.kfold):
                datasets = {
                    'training': (self.x_train[fold], self.y_train[fold]),
                    'validation': (self.x_val[fold], self.y_val[fold]),
                    'test': (self.x_test, self.y_test)
                }
                datasets.update(self.add_test if self.add_test else {})

                for ds_name, (x, y) in datasets.items():
                    if ds_name not in ('training', 'validation', 'test'):
                        x = self.scaler.transform(x)
                    start = perf_counter()
                    prediction = models[fold].predict(x)
                    time = perf_counter() - start
                    n_x = len(x)
                    rel_time = time / n_x

                    metric_arr = calc_metrics(y.ravel(), prediction.ravel())
                    q2_f3 = q2_f3_score(y, prediction, self.y_train[fold])
                    self._predictions2frame(model_name, ds_name, prediction, metric_arr, q2_f3, rel_time)
                    if self.n_perm:
                        for feat_idx in range(self.n_features):
                            for perm in range(self.n_perm):
                                x_cop = x.copy()
                                np.random.shuffle(x_cop[:, feat_idx])
                                y_perm = models[fold].predict(x_cop)
                                self._fi2frame(model_name, ds_name, fold, feat_idx, perm, y, y_perm, metric_arr, q2_f3)

            for ds_name in self.score_dict:
                self.score_dict[ds_name][model_name] = pd.DataFrame.from_dict(self.score_dict[ds_name][model_name])
                self.score_dict[ds_name][model_name].index.name = "fold"
        self.score_dict = {ds_name: pd.concat(dat, names=["model"]) for ds_name, dat in self.score_dict.items()}
        out = pd.concat(self.score_dict, names=["dataset"])
        if self.n_perm:
            mindex = pd.MultiIndex.from_tuples(self.perm_metrics.keys(),
                                               names=('model', 'dataset', 'fold', 'feature', 'permutation'))
            perm_df = pd.DataFrame([*self.perm_metrics.values()], index=mindex,
                                   columns=('mean_squared_error',
                                            'median_absolute_error',
                                            'r_spearman',
                                            'r_pearson',
                                            'r2_score',
                                            'Q2_F3'))
            return out, perm_df
        return out

    def _predictions2frame(self, model_name, dataset, pred_arr, metric_arr, q2_f3, rel_time):
        mse, mae, rho, r, r2 = metric_arr
        score_dict = self.score_dict.setdefault(dataset, {})
        score_dict.setdefault(model_name, {})
        score_dict[model_name].setdefault('mean_squared_error', []).append(mse)
        score_dict[model_name].setdefault('median_absolute_error', []).append(mae)
        score_dict[model_name].setdefault("RMSE", []).append(np.sqrt(mse))
        score_dict[model_name].setdefault('r_pearson', []).append(r)
        score_dict[model_name].setdefault('r_spearman', []).append(rho)
        score_dict[model_name].setdefault('r2_score', []).append(r2)
        score_dict[model_name].setdefault('Q2_F3', []).append(q2_f3)
        score_dict[model_name].setdefault('raw_predictions', []).append(pred_arr)
        score_dict[model_name].setdefault('time', []).append(rel_time)

    def _fi2frame(self, model_name, dataset, fold, feature_idx, perm_idx,
                  y_true, y_perm, baseline_metrics, baseline_q2_f3):
        alt_metrics = calc_metrics(y_true.ravel(), y_perm.ravel())
        alt_q2_f3 = q2_f3_score(y_true, y_perm, self.y_train[fold])
        delta_metrics = [alt-base for alt, base in zip(alt_metrics, baseline_metrics)]
        delta_q2_f3 = alt_q2_f3 - baseline_q2_f3
        delta_metrics.append(delta_q2_f3)
        if self.feature_names:
            feature_idx = self.feature_names[feature_idx]
        self.perm_metrics[(model_name, dataset, fold, feature_idx, perm_idx)] = delta_metrics


def calc_metrics(true, predicted):
    """
    :param true:
    :param predicted:
    :return: mse, mae, rho, r, r2
    """
    mse = metrics.mean_squared_error(true, predicted)
    mae = metrics.median_absolute_error(true, predicted)
    r2 = metrics.r2_score(true, predicted)
    rho, _ = stats.spearmanr(true, predicted)
    r, _ = stats.pearsonr(true, predicted)
    return mse, mae, rho, r, r2


def find_best_ml_models(candidate_dict, x_train, y_train, x_val, y_val):
    """
    Hyperparameter search for sklearn models
    :param candidate_dict:
    :param x_train:
    :param y_train:
    :param x_val:
    :param y_val:
    :return: trained_models, hyper_params (dicts)
    """
    trained_models = {}
    hyper_params = {}
    for candidate_name, candidate in candidate_dict.items():
        regr_class, regr_params = candidate
        regr_param_list = [*ParameterGrid(regr_params)]

        def score(model, x, y):  # Place the appropriate metric here
            return stats.pearsonr(y, model.predict(x))[0]

        def eval_ml_model(params):
            models = [regr_class(**params).fit(x, y) for x, y in zip(x_train, y_train)]
            return np.mean([score(model, x, y) for model, x, y in zip(models, x_val, y_val)]), models

        scores, _models = zip(*[eval_ml_model(par) for par in regr_param_list])

        best_index = int(np.argmax(scores))
        hyper_params[candidate_name] = (regr_class, regr_param_list[best_index])
        trained_models[candidate_name] = _models[best_index]

    return trained_models, hyper_params
