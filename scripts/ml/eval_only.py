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

# Utility to run analysis including feature importance on existing model trainings and additional external datasets
# Author(s): Stefan Holderbach
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import os
import warnings
import shutil
import argparse
import pandas as pd
import numpy as np
from sklearn.externals import joblib
from keras.models import load_model
from src.manage import ManagementSystem
from src.loader import DataLoader
from src.MLA import MLA

parser = argparse.ArgumentParser(
    description='Predict on datasets, compute metrics and if necessary feature importances')
parser.add_argument("-i", "--indir", type=str, required=True,
                    help="please set directory with the results from a training run")
parser.add_argument("-o", "--outdir", type=str, required=True,
                    help="please set directory. where data from evaluation will be written")
parser.add_argument("-e", "--testfiles", type=str, required=False, nargs='*',
                    help="Additional files for testing")
parser.add_argument("-fi", "--feature_importance", action="store_true", default=False,
                    help="Set to calculate feature importance")

args = parser.parse_args()

## PLACE your model names from train_and_eval.py here
ml_names = [
    "Mean of training set",
    "Nearest-Neighbor",
    "Support Vector Regression",
    "Epsilon-Support Vector Regression",
    "LinearRegression",
    "Random Forest",
    "Extremely Random Forest"]
dl_names = ["RASPDeep"]

emergency_rename = {"PA (D+E)": "PA(D+E)", "PMR (Arom)": "PMR(Arom)"}

def ml_parallelize(model):
    if "n_jobs" in model.get_params():
        model.set_params(n_jobs=-1)
    return model


if __name__ == '__main__':

    # CONST
    ds_modes = ("val", "test", "train")
    # Load configuration
    outdir = args.outdir
    indir = args.indir
    ext_test_files = args.testfiles

    management_system = ManagementSystem(datadir=outdir, make_unique=False)
    wdir, logger = management_system.setting()

    logger.info(f"Training taken from {indir}")
    logger.info(f"External test files: {args.testfiles}")
    logger.info(f"Calculate feature importance: {args.feature_importance}")

    # Load the training related data, scalers and used features
    (index_train, x_train, y_train,
     index_val, x_val, y_val,
     index_test, x_test, y_test) = [
        np.load(os.path.join(indir, f'{fname}.npy'), allow_pickle=True) for fname in
        ('index_train', 'x_train', 'y_train',
         'index_val', 'x_val', 'y_val',
         'index_test', 'x_test', 'y_test')
    ]

    scalers = joblib.load(os.path.join(indir, 'scalers.pkl'))
    feature_names = joblib.load(os.path.join(indir, 'feature_names.pkl'))
    feature_names = [emergency_rename[name] if name in emergency_rename else name 
                      for name in feature_names
                    ] 

    # Drop the included atom efficiency
    y_train, y_val, y_test = [yarr.take(0, axis=-1) for yarr in (y_train, y_val, y_test)]

    # Extract further necessary information

    n_shuffle = x_train.shape[0]
    n_folds = x_train.shape[1]

    logger.info(f"Training had been performed with {n_shuffle} replicates and {n_folds}-fold cv")
    logger.info(f"Features considered: {feature_names}")

    # Load external test data
    ext_test_data = None
    if ext_test_files:
        ext_test_data = {}
        ext_test_files = DataLoader(indir, ext_test_files, logger,
                                    to_shuffle=False, ).retrieve_data()
        for filename, dataset in ext_test_files.items():
            try:
                ds_name = "_".join(os.path.basename(filename).split(".")[:-1])
                x = dataset[feature_names]
                y = dataset["Expt_BE"]
                ext_test_data[ds_name] = (x, y)
            except:
                logger.exception(f'Problem with {filename}')

    # Prepare data structures for the results here
    metrics_dfs = {}
    fi_dfs = {}

    # The big loop
    for i in range(n_shuffle):
        # Load the models
        model_dict = {name: [ml_parallelize(joblib.load(os.path.join(indir, f"{name}_fold_{j}_shuffle_{i}.joblib")))
                       for j in range(n_folds)]
                      for name in ml_names}
        model_dict.update({
                        name: [load_model(os.path.join(indir, f"{name}_fold_{j}_shuffle_{i}.h5"))
                       for j in range(n_folds)]
                      for name in dl_names
        })


        kwargs = {}
        if args.feature_importance:
            kwargs = {"n_perm": 5,
                      "feature_names": feature_names}

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
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
    # Make raw csv output for humans or bash use
    os.makedirs(os.path.join(wdir, "raw_output"), exist_ok=True)
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

    # Dump the rest of useful information
    joblib.dump(metric_data.raw_predictions, os.path.join(wdir, "raw_predictions.joblib"))
    metric_data = metric_data.drop(columns="raw_predictions")
    metric_data.to_csv(os.path.join(wdir, "all_metrics.tsv"), sep="\t")
    joblib.dump(metric_data, os.path.join(wdir, "all_metrics.joblib"))
    avg_metrics = metric_data.groupby(level=["model", "dataset"]).mean()
    std_metrics = metric_data.groupby(level=["model", "dataset"]).std()
    avg_metrics.to_csv(os.path.join(wdir, "avg_metrics.tsv"), sep="\t")
    std_metrics.to_csv(os.path.join(wdir, "std_metrics.tsv"), sep="\t")

    if args.feature_importance:
        fi_data = pd.concat(fi_dfs, names=["replicate"])
        joblib.dump(fi_data, os.path.join(wdir, "feature_importances.joblib"))
        fi_data.to_csv(os.path.join(wdir, "feature_importances.tsv"), sep="\t")

    # As a convenience for the existing plotting pipeline copy the hyperparameter info to the new directory
    try:
        shutil.copy(os.path.join(indir, "hyper.tsv"), os.path.join(wdir, "hyper.tsv"))
    except FileNotFoundError:
        pass
    with open(os.path.join(wdir, "my_training.txt"), "w") as txt_file:
        txt_file.write(os.path.abspath(indir))
