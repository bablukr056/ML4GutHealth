#! /lustrehome/boss/miniconda3/envs/ML_models/bin/python

import os
import sys
import argparse
import argcomplete



def be_parser():
    parser = argparse.ArgumentParser(
        description="evaluate SVC models",
        prefix_chars="-")
    parser.add_argument("-m", "--model", help="model path",
                        action="store", required=True, type=str)
    argcomplete.autocomplete(parser)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def shap_on_models(model_path):
    import joblib
    import shap
    import pandas as pd
    model = joblib.load(model_path)
    # importiamo il dataste di training
    X_train = pd.read_csv("train_dataset.csv", header=0, index_col=0, sep=",")
    Y_train = X_train["Class_Label"].values
    X_train.drop(["Class_Label"], axis=1, inplace=True)
    # importiamo il dataset di test
    X_test = pd.read_csv("test_dataset.csv", header=0, index_col=0, sep=",")
    Y_test = X_test["Class_Label"]
    X_test.drop(["Class_Label"], axis=1, inplace=True)
    background = shap.kmeans(X_train, 100)
    explainer = shap.KernelExplainer(
        model.predict_proba,
        background,
        algorithm="permutation"
    )
    import numpy as np

    X_explain = X_test.sample(
        n=500,
        random_state=42
    )

    chunks = np.array_split(X_explain, 100)

    results = joblib.Parallel(
        n_jobs=100,
        backend="loky"
    )(
        joblib.delayed(explainer.shap_values)(
            chunk,
            nsamples=50
        )
        for chunk in chunks
    )
    shap_values = np.concatenate(results, axis=0)
    outfile = model_path.replace(".joblib", ".shap.joblib")
    joblib.dump(shap_values, outfile)

if __name__ == "__main__":
    arg = be_parser()
    model_path  = arg.model
    # with open(f"{rank}_prevalence_exec.log", "w") as log:
    print(arg.__dict__)
    shap_on_models(model_path)
