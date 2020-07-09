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

# RASPD+ inference script to run predictions from a saved models folder
# Author(s): Stefan Holderbach
# Licensed under the EUPL 1.2
# Code for "RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features"
import os
import numpy as np
import argparse
import logging
import tensorflow as tf
from sklearn.externals import joblib
from keras.models import load_model
from keras import backend as keras_backend
from src.loader import DataLoader
from multiprocessing import cpu_count


possible_models = {'esvr': 'Epsilon-Support Vector Regression',
                   'erf': 'Extremely Random Forest',
                   'lr': 'LinearRegression',
                   'knn': 'Nearest-Neighbor',
                   'dnn': 'RASPDeep',
                   'rf': 'Random Forest',
                   'svr': 'Support Vector Regression'}

emergency_rename = {"PA (D+E)": "PA(D+E)", "PMR (Arom)": "PMR(Arom)"}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Predict binding energy from protein-ligand complex features')
    parser.add_argument("-d", "--datadir", type=str, required=True, help="please set output directory")
    parser.add_argument("-w", "--wdir", type=str, required=True, help="please set weights directory")
    parser.add_argument("-i", "--inputfile", type=str, required=True,
                        help="please specify input file (either excel, csv or text)")
    parser.add_argument("-m", "--models", type=str, nargs='*', default=["erf"],
                        help=f"please specify at least one model from {possible_models}")
    parser.add_argument("-c", "--cpu_count", type=int, default=None, help="Number of CPU cores to use")
    parser.add_argument("-s", "--select_shuffle", type=int, default=0, help="Which of the shuffles should be used?")
    args = parser.parse_args()

    datadir, wdir, filename, models = args.datadir, args.wdir, args.inputfile, args.models
    uninferred = True

    os.makedirs(datadir, exist_ok=True)

    # Logging setup
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=os.path.join(datadir, 'logfile.log'), level=logging.INFO)

    # Limit core count for keras
    logger.info(f"Core count limit={args.cpu_count}")
    if args.cpu_count:
        if args.cpu_count > cpu_count():
            args.cpu_count = cpu_count()
        tf_config = tf.ConfigProto(intra_op_parallelism_threads=args.cpu_count,
                                   inter_op_parallelism_threads=args.cpu_count,
                                   allow_soft_placement=True,
                                   device_count={'CPU': args.cpu_count})
        session = tf.Session(config=tf_config)
        keras_backend.set_session(session)

    file = os.path.join(datadir, filename)

    try:
        y_train = np.load(os.path.join(wdir, 'y_train.npy'))
        feature_names = joblib.load(os.path.join(wdir, 'feature_names.pkl'))
        feature_names = [emergency_rename[name] if name in emergency_rename else name 
                          for name in feature_names
                        ] 

        assert y_train.ndim == 4, "Stored trainings data set does not seem like nested crossvalidation"
        n_folds = y_train.shape[1]  # HERE we define the number of expected k-folds that will be averaged
        shuffle = args.select_shuffle  # Which shuffle to take
        assert 0 <= shuffle < y_train.shape[0], "-s/--select_shuffle option not within range of trained replicates"

        scaler = joblib.load(os.path.join(wdir, 'scalers.pkl'))[shuffle]

    except FileNotFoundError:
        logger.error("Stored objects not found; are you sure this is a weight directory?")
        exit(1)
    except AssertionError as err:
        logger.error(err.args)
        exit(1)

    try:
        dl = DataLoader(datadir=datadir, filenames=[filename], logger=logger, )
        data = dl.retrieve_data()[filename]
        data = data[feature_names]
        x = scaler.transform(data)
    except:
        logger.error("An error occured trying to load the input datafile. "
                     "Please check the file location and format according to the spec.")
        logger.exception("ABORTING")
        exit(1)

    for model_name in models:
        preds = []
        for fold in range(n_folds):
            mname = possible_models.get(model_name)
            basename = os.path.join(wdir, f'{mname}_fold_{fold}_shuffle_{shuffle}')
            if model_name == "dnn":
                model = load_model(basename + '.h5')
            elif model_name in possible_models:
                model = joblib.load(basename + '.joblib')
                if "n_jobs" in model.get_params():
                    model.set_params(n_jobs=(args.cpu_count or -1))
            else:
                logger.warning(f"Specified model name {model_name} not found")
                break
            preds.append(model.predict(x))

        else:  # Only executed if folds have been properly read
            preds = np.reshape(preds, (n_folds, -1))
            pred_be = np.stack([np.mean(preds, axis=0), np.std(preds, axis=0)], axis=1)
            np.savetxt(os.path.join(datadir, f'{model_name}.out'), pred_be, fmt="%.4f\t%.4f")
            logger.info(f"Succesfully inferred with {model_name}")
            logger.info(f'Written results to {os.path.join(datadir, model_name)}.out')
            uninferred = False

    if uninferred:
        logger.warning("No inference on any model took place!")
        exit(1)
    logger.info("FINISHED")
