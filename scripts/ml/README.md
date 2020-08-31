# Machine learning submodule

This folder contains the code to train end evaluate the RASPD+ machine learning models

## Content

* `/src`
  * python modules that drive the model training and load and prepare the datasets
* `/notebooks`
  * `plot_results.ipynb`
    * makes all the figures for the analysis of the machine learning models seen in the paper and also provides additional information
  * `datastet_stats.ipynb`
    * analyzes the properties of our chosen descriptors on the PDBbind data set as well as the external test sets.
    * Also provides figures found in the paper
  * `*_external_benchmarks.csv` 
    * Benchmark results from other methods as seen in Jiménez et al. 2018
  * `verify_old_RASPD.ipynb`
    * read results data from the old RASPD approach (Mukherjee et al. 2011)
  * `dud_e_enrichment_assessment.ipynb`
    * Compute enrichment factors for evaluation on the Directory of useful Decoys - Enhanced (DUD-E) data set
  * `launch_runs.sh`
    * How to run the `train_and_eval.py` script to obtain the models used in `plot_results.ipynb`
* `pdbbind_descriptors.xlsx`
  * Precomputed RASPD+ descriptor values from the PDBbind refined set
* `dataset_smaller.txt`
  * small dataset from the PDBbind set for debugging
* `train_and_eval.py`
  * script to train new models and directly perform evaluation on selected datasets.
  * Modify source code to include additional hyperparameter options
* `eval_only.py`
  * Takes the result folder from a `train_and_eval.py` run and performs evaluation or feature importance analysis on additional datasets
  * output folders are consumed by the `plot_results` notebook
* `infer.py`
  * Script to predict binding free energies for a set of computed descriptors with an existing model.
  * Performs bagging by combining the results from cross-validation folds into one prediction
  
For information on how to run the three command line scripts, run them with the `--help` option.

## `train_and_eval.py`

Train classifier to predict binding energy from protein-ligand complex
features

```
arguments:
  -h, --help            show this help message and exit
  -d DATADIR, --datadir DATADIR
                        please set data directory. Where to write data (if -p
                        not set creates a unique dir inside)
  -f [TRAINFILES [TRAINFILES ...]], --trainfiles [TRAINFILES [TRAINFILES ...]]
                        please set file for training
  -e [TESTFILES [TESTFILES ...]], --testfiles [TESTFILES [TESTFILES ...]]
                        Additional files for testing
  -c [COLUMNS_TO_IGNORE [COLUMNS_TO_IGNORE ...]], --columns_to_ignore [COLUMNS_TO_IGNORE [COLUMNS_TO_IGNORE ...]]
                        (optional) ignore feature columns
  -t TEST_SIZE, --test_size TEST_SIZE
                        optional: set test size
  -j OUTER_CV_SPLITS, --outer_cv_splits OUTER_CV_SPLITS
                        optional: set random test set splits
  -k INNER_CV_SPLITS, --inner_cv_splits INNER_CV_SPLITS
                        optional: set inner CV kfold splits
  -l LOADING, --loading LOADING
                        please specify if dataset shall be loaded from a
                        previous run
  -fi, --feature_importance
                        Set to calculate feature importance
  -p, --path_fixed      Set to use the datadir directly and don't use a unique
                        directory
```

## `eval_only.py`

Predict on datasets, compute metrics and if necessary feature importances

```
arguments:
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        please set directory with the results from a training
                        run
  -o OUTDIR, --outdir OUTDIR
                        please set directory. where data from evaluation will
                        be written
  -e [TESTFILES [TESTFILES ...]], --testfiles [TESTFILES [TESTFILES ...]]
                        Additional files for testing
  -fi, --feature_importance
                        Set to calculate feature importance
```

## `infer.py`

Predict binding energy from protein-ligand complex features

```
arguments:
  -h, --help            show this help message and exit
  -d DATADIR, --datadir DATADIR
                        please set output directory
  -w WDIR, --wdir WDIR  please set weights directory
  -i INPUTFILE, --inputfile INPUTFILE
                        please specify input files (either excel, csv or text)
  -m [MODELS [MODELS ...]], --models [MODELS [MODELS ...]]
                        please specify at least one model from {'esvr':
                        'Epsilon-Support Vector Regression', 'erf': 'Extremely
                        Random Forest', 'lr': 'LinearRegression', 'knn':
                        'Nearest-Neighbor', 'dnn': 'RASPDeep', 'rf': 'Random
                        Forest', 'svr': 'Support Vector Regression'}
  -c CPU_COUNT, --cpu_count CPU_COUNT
                        Number of CPU cores to use
  -s SELECT_SHUFFLE, --select_shuffle SELECT_SHUFFLE
                        Which of the shuffles should be used?
```

## License

All code in this section is licensed under the European Union Public License v. 1.2
