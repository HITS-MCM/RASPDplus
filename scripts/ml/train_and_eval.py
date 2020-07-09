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

# RASPD+ training and evaluation script
# Author(s): Stefan Holderbach, Lukas Adam
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import os
import argparse
import pandas as pd
import numpy as np
from sklearn.externals import joblib
from src.manage import ManagementSystem
from src.prep import DataPrep
from src.loader import DataLoader
from src.model import DL, BestModelKCross
from src.MLA import MLA, find_best_ml_models
from sklearn import linear_model, neighbors, svm, ensemble, dummy

parser = argparse.ArgumentParser(
    description='Train classifier to predict binding energy from protein-ligand complex features')
parser.add_argument("-d", "--datadir", type=str, required=True,
                    help="please set data directory. Where to write data (if -p not set creates a unique dir inside)")
parser.add_argument("-f", "--trainfiles", type=str, required=True, nargs='*',
                    help="please set file for training")
parser.add_argument("-e", "--testfiles", type=str, required=False, nargs='*',
                    help="Additional files for testing")
parser.add_argument("-c", "--columns_to_ignore", type=str, required=False, nargs='*',
                    help="(optional) ignore feature columns")
parser.add_argument("-t", "--test_size", type=float,
                    required=False, default=0.125, help="optional: set test size")
parser.add_argument("-j", "--outer_cv_splits", type=int,
                    required=False, default=10, help="optional: set random test set splits")
parser.add_argument("-k", "--inner_cv_splits", type=int,
                    required=False, default=6, help="optional: set inner CV kfold splits")
parser.add_argument("-l", "--loading", type=bool, default=False, required=False,
                    help="please specify if dataset shall be loaded from a previous run")
parser.add_argument("-fi", "--feature_importance", action="store_true", default=False,
                    help="Set to calculate feature importance")
parser.add_argument("-p", "--path_fixed", action="store_true", default=False,
                    help="Set to use the datadir directly and don't use a unique directory")

args = parser.parse_args()

if __name__ == '__main__':

    # CONST
    ds_modes = ("val", "test", "train")
    # Load configuration
    nn_param_space = {
        "lr": [0.1, 0.3, 1],  # Attention: For true SGD learning rate is multiplied by 0.01
        'batch_size': [8],
        'epochs': [250],
        'dropout1': [0.1],
        'dropout2': [0.1],
        'loss_fn': ["mean_squared_error"],
        'opt': ['sgd'],
        'activation1': ["elu"],
        'activation2': ["elu"],
        'node1': [20],
        'node2': [10],
        'initializer': ["VarianceScaling", "RandomNormal"]
    }
    ml_regressor_classes = {
        "Mean of training set":
            (dummy.DummyRegressor, {}),
        "Nearest-Neighbor":
            (neighbors.KNeighborsRegressor, {"n_neighbors": [1, 5, 10, 20]}),
        "Support Vector Regression":
            (svm.LinearSVR, {"C": [0.01, 0.1, 1, 10],
                             "epsilon": [0.0]}),
        "Epsilon-Support Vector Regression":
            (svm.SVR, {"C": [0.01, 0.1, 1, 10],
                       "epsilon": [0.1],
                       "kernel": ["rbf"]}),
        "LinearRegression":
            (linear_model.LinearRegression, {}),
        "Random Forest":
            (ensemble.RandomForestRegressor, {"n_estimators": [20, 100, 200],
                                              "max_features": ["log2", "sqrt"],
                                              "min_samples_leaf": [1, 3, 5]}),
        "Extremely Random Forest":
        # https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.ExtraTreesRegressor.html
            (ensemble.ExtraTreesRegressor, {"n_estimators": [20, 100, 200],
                                            "max_features": ["log2", "sqrt"],
                                            "min_samples_leaf": [1, 3, 5]})
    }

    datadir = args.datadir
    trainfiles = args.trainfiles
    ext_test_files = args.testfiles
    n_splits = args.inner_cv_splits
    n_shuffle = args.outer_cv_splits
    test_size = args.test_size
    ignore_cols = args.columns_to_ignore

    management_system = ManagementSystem(datadir=datadir, make_unique=(not args.path_fixed))
    wdir, logger = management_system.setting()

    logger.info(f"Specified random test ratio: {args.test_size}")
    logger.info(f"Number of random test splits: {args.outer_cv_splits}")
    logger.info(f"Number of inner cross-validation splits: {args.inner_cv_splits}")
    logger.info(f"File for train/test/val: {args.trainfiles}")
    logger.info(f"External test files: {args.testfiles}")
    logger.info(f"Calculate feature importance: {args.feature_importance}")
    logger.info(f"The following features will be ignored: {ignore_cols}")

    # Load external test data
    ext_test_data = None
    if ext_test_files:
        ext_test_data = {}
        ext_test_files = DataLoader(datadir, ext_test_files, logger,
                                    to_shuffle=False, ignore_cols=ignore_cols).retrieve_data()
        for filename, dataset in ext_test_files.items():
            ds_name = "_".join(os.path.basename(filename).split(".")[:-1])
            x = dataset.drop(['Expt_BE', 'Atom Efficiency', 'No. of atoms'], axis=1)
            y = dataset["Expt_BE"]
            ext_test_data[ds_name] = (x, y)

    # Load data either from files and perform splitting or use prepared datasets
    if args.loading is False:
        data_loader = DataLoader(datadir, trainfiles, logger, to_shuffle=False, ignore_cols=ignore_cols)
        data = data_loader.retrieve_data()
        data = data[trainfiles[0]]

        # extract expected binding energies from the dataset
        be_and_ae = data[['Expt_BE', "Atom Efficiency"]]
        data = data.drop(['Expt_BE', 'Atom Efficiency', 'No. of atoms'], axis=1, errors="ignore")

        data.columns = [f.strip() for f in list(data.columns)]
        feature_names = data.columns
        print('\ndata shape:', data.shape)

        indices = np.array(data.index)
        x = np.array(data)
        y = np.array(be_and_ae)  # Atom efficiency is bundled in for filtering and analysis (will be removed)

        data_prep = DataPrep((x, y, indices), test_size=test_size, n_splits=n_splits,
                             n_shuffle=n_shuffle, single_mode=False)
        data_prep.save_datasets(wdir)
        train_val_test, scalers = data_prep.retrieve_datasets()
        joblib.dump(feature_names, os.path.join(wdir, "feature_names.pkl"))
        joblib.dump(scalers, os.path.join(wdir, 'scalers.pkl'))

        (index_train, x_train, y_train,
         index_val, x_val, y_val,
         index_test, x_test, y_test) = [
            train_val_test[fname] for fname in
            ('index_train', 'x_train', 'y_train',
             'index_val', 'x_val', 'y_val',
             'index_test', 'x_test', 'y_test')]
    else:
        (index_train, x_train, y_train,
         index_val, x_val, y_val,
         index_test, x_test, y_test) = [
            np.load(os.path.join(indir, f'{fname}.npy'), allow_pickle=True) for fname in
            ('index_train', 'x_train', 'y_train',
             'index_val', 'x_val', 'y_val',
             'index_test', 'x_test', 'y_test')
        ]
        scalers = joblib.load(os.path.join(datadir, 'scalers.pkl'))
        feature_names = joblib.load(os.path.join(datadir, 'feature_names.pkl'))

    y_train, y_val, y_test = [yarr.take(0, axis=-1) for yarr in (y_train, y_val, y_test)]

    # Prepare data structures for the results here
    metrics_dfs = {}
    fi_dfs = {}
    best_hyper = []

    # The big loop
    for i in range(n_shuffle):
        # Find appropriate hyperparameters DL
        dl = DL(x_train[i], y_train[i], x_val[i], y_val[i], nn_param_space, wdir)
        taloslib = dl.run()
        talosruns = [scan['talosScan'] for scan in taloslib.values()]
        dl_kcross = BestModelKCross(talosruns)
        best_dl_model = dl_kcross.get_saved_models()
        best_dl_params = dl_kcross.params
        rundir = 'deep_learning'
        for model_name, model in best_dl_model.items():
            model.save(os.path.join(wdir, f'{model_name}_shuffle_{i}.h5'))

        # Find hyperparameters ML
        ml_models, ml_best_params = find_best_ml_models(ml_regressor_classes,
                                                        x_train[i], y_train[i],
                                                        x_val[i], y_val[i])
        for name, models in ml_models.items():
            for j, model in enumerate(models):
                fname = os.path.join(wdir, f"{name}_fold_{j}_shuffle_{i}.joblib")
                joblib.dump(model, fname)
        ml_best_params["RASPDeep"] = best_dl_params
        best_hyper.append(ml_best_params)

        # Evaluate DL and ML
        # Evaluate the following things: Everything on the normal train/val/test, external test sets
        # Run on the column wise shuffled test sets
        model_dict = ml_models
        model_dict["RASPDeep"] = [*best_dl_model.values()]

        kwargs = {}
        if args.feature_importance:
            kwargs = {"n_perm": 5,
                      "feature_names": [*feature_names]}

        mla = MLA(x_train[i], y_train[i],
                  x_val[i], y_val[i],
                  x_test[i], y_test[i],
                  additional_tests=ext_test_data, scaler=scalers[i], **kwargs)
        if not args.feature_importance:
            ml_metrics = mla.run(model_dict)
        else:
            ml_metrics, fi_df = mla.run(model_dict)
            fi_dfs[i] = fi_df
        metrics_dfs[i] = ml_metrics

    metric_data = pd.concat(metrics_dfs, names=["replicate"])
    os.makedirs(os.path.join(wdir, "raw_output"), exist_ok=True)
    # TODO: insert the human readable prediction export here
    for index_tup, preds in metric_data["raw_predictions"].items():
        repl, ds, mdl_name, fold = index_tup
        index = None
        trues = None
        preds = np.ravel(preds)
        if ds == "test":
            index = index_test[repl]
            trues = y_test[repl]
        elif ds == "validation":
            index = index_val[repl, fold]
            true = y_val[repl, fold]
        elif ds == "training":
            index = index_val[repl, fold]
            true = y_val[repl, fold]
        else:
            index = ext_test_data[ds][0].index
            trues = ext_test_data[ds][1].values
        out_fname = f"{ds}_{mdl_name}_{repl}_{fold}.csv".replace(" ", "_")
        out_df = pd.DataFrame({"true": trues, "pred": preds})
        out_df.reindex(index)
        out_df.to_csv(os.path.join(wdir, "raw_output", out_fname), sep=";", index_label="PDBID")

    joblib.dump(metric_data.raw_predictions, os.path.join(wdir, "raw_predictions.joblib"))
    metric_data = metric_data.drop(columns="raw_predictions")
    metric_data.to_csv(os.path.join(wdir, "all_metrics.tsv"), sep="\t")
    joblib.dump(metric_data, os.path.join(wdir, "all_metrics.joblib"))
    avg_metrics = metric_data.groupby(level=["model", "dataset"]).mean()
    std_metrics = metric_data.groupby(level=["model", "dataset"]).std()
    avg_metrics.to_csv(os.path.join(wdir, "avg_metrics.tsv"), sep="\t")
    std_metrics.to_csv(os.path.join(wdir, "std_metrics.tsv"), sep="\t")

    print(pd.DataFrame(best_hyper))
    pd.DataFrame(best_hyper).to_csv(os.path.join(wdir, "hyper.tsv"), sep="\t", index=False)

    if args.feature_importance:
        fi_data = pd.concat(fi_dfs, names=["replicate"])
        joblib.dump(fi_data, os.path.join(wdir, "feature_importances.joblib"))
        fi_data.to_csv(os.path.join(wdir, "feature_importances.tsv"), sep="\t")
