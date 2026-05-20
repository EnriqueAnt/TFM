#!/usr/bin/env python
# coding: utf-8

import os

# 1. CONFIGURACIÓN DE ENTORNO (Debe ir antes de importar numpy/sklearn)
# Bloquea el paralelismo desde la CPU para evitar variaciones decimales en los cálculos
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import json
import random
import numpy as np
import pandas as pd

# 2. FUNCIÓN DE SEMILLA GLOBAL
def set_reproducibility(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    # Si en el futuro usas bibliotecas como PyTorch o TensorFlow, 
    # sus semillas se fijarían aquí también.

# Ejecutar fijación de semilla al inicio
set_reproducibility(42)

# 3. IMPORTACIONES RESTANTES
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    balanced_accuracy_score,
    precision_score,
    f1_score,
    accuracy_score
)
from sklearn.cluster import MiniBatchKMeans, KMeans

from imblearn.under_sampling import (
    ClusterCentroids, CondensedNearestNeighbour, EditedNearestNeighbours,
    RepeatedEditedNearestNeighbours, AllKNN, NearMiss,
    NeighbourhoodCleaningRule, OneSidedSelection,
    RandomUnderSampler, TomekLinks
)
from imblearn.over_sampling import (
    RandomOverSampler, SMOTE,
    BorderlineSMOTE, KMeansSMOTE
)
from imblearn.combine import SMOTEENN, SMOTETomek

# --- BLOQUE DE ARGUMENTOS ---
parser = argparse.ArgumentParser(
    description='Experimento 1.3 (R3): Resampling -> Feature Selection'
)
parser.add_argument('--dataset', type=str, required=True)
parser.add_argument('--k', type=int, required=True)
parser.add_argument('--input_dir', type=str, default='.')
parser.add_argument('--output_root', type=str, default='results_ml_classification')

args = parser.parse_args()

# --- FIX: Asegurar rutas absolutas al directorio Project ---
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_dir = os.path.dirname(_script_dir)
if getattr(args, 'output_root', None) == 'results_ml_classification':
    args.output_root = os.path.join(_project_dir, 'results_ml_classification')
if getattr(args, 'input_dir', None) == '.':
    args.input_dir = _project_dir
elif getattr(args, 'input_dir', None) == 'results_ml_classification':
    args.input_dir = os.path.join(_project_dir, 'results_ml_classification')


# --- CONFIGURACIÓN DINÁMICA ---
file_name = f"genes-{args.dataset}"
n_best_features = args.k
path_workflow = "R3"
random_state_base = 42
n_runs = 10

# Carpetas
current_output_dir = os.path.join(
    args.output_root,
    args.dataset,
    f"k_{n_best_features}"
)
os.makedirs(current_output_dir, exist_ok=True)

metadata_dir = os.path.join(current_output_dir, "runs_metadata")
os.makedirs(metadata_dir, exist_ok=True)

print(f"--- INICIO SCRIPT 1.3 (R3) ---")
print(f"Dataset: {file_name} | K: {n_best_features}")

# --- CARGA DE DATOS ---
input_path = os.path.join(args.input_dir, f"{file_name}.csv")
df = pd.read_csv(input_path)

# Split
X = pd.get_dummies(df.iloc[:, :-1])
y = df.iloc[:, -1]

X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    stratify=y,
    random_state=random_state_base
)

# --- BASELINE (sin resampling) ---
selector_init = SelectKBest(f_classif, k=n_best_features)
X_train_selected_init = selector_init.fit_transform(X_train, y_train)
X_test_selected_init = selector_init.transform(X_test)

selected_features_init = X_train.columns[selector_init.get_support()]

# Guardado referencia
train_df = pd.DataFrame(X_train_selected_init, columns=selected_features_init)
train_df["target"] = y_train.values
test_df = pd.DataFrame(X_test_selected_init, columns=selected_features_init)
test_df["target"] = y_test.values

train_df.to_csv(
    os.path.join(current_output_dir,
    f"{file_name}-resampling_FR-{n_best_features}_train.csv"),
    index=False
)

test_df.to_csv(
    os.path.join(current_output_dir,
    f"{file_name}-resampling_FR-{n_best_features}_test.csv"),
    index=False
)

# --- PARÁMETROS DINÁMICOS ---
if args.dataset == "cogdx":
    k_neighbors_val = 2
    n_clusters_val = 1
    balance_threshold_val = 0.001
else:
    k_neighbors_val = 3
    n_clusters_val = 2
    balance_threshold_val = 0.01

# --- MÉTODOS DE RESAMPLING ---
base_methods = {
    "CC": ClusterCentroids(
        estimator=MiniBatchKMeans(n_init="auto", random_state=random_state_base)
    ),
    "CNN": CondensedNearestNeighbour(),
    "ENN": EditedNearestNeighbours(),
    "RENN": RepeatedEditedNearestNeighbours(),
    "AllKNN": AllKNN(),
    "NearMiss": NearMiss(),
    "NCR": NeighbourhoodCleaningRule(),
    "OSS": OneSidedSelection(),
    "RUS": RandomUnderSampler(),
    "Tomek": TomekLinks(),
    "ROS": RandomOverSampler(),
    "SMOTE": SMOTE(k_neighbors=k_neighbors_val),
    "BorderlineSMOTE": BorderlineSMOTE(k_neighbors=3, m_neighbors=3),
    "KMeansSMOTE": KMeansSMOTE(
        cluster_balance_threshold=balance_threshold_val,
        kmeans_estimator=KMeans(
            n_init="auto",
            n_clusters=n_clusters_val,
            random_state=random_state_base
        )
    ),
    "SMOTEENN": SMOTEENN(smote=SMOTE(k_neighbors=k_neighbors_val)),
    "SMOTETomek": SMOTETomek(smote=SMOTE(k_neighbors=k_neighbors_val)),
}

# RF independiente por método
resampling_methods = {
    f"{name}_RF": (
        sampler,
        RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            min_samples_leaf=5,
            min_samples_split=10,
            max_features="sqrt",
            class_weight=None
        )
    )
    for name, sampler in base_methods.items()
}

# --- ESTRUCTURA DE MÉTRICAS ---
metrics = {
    "BA": balanced_accuracy_score,
    "PS": lambda y_true, y_pred:
        precision_score(y_true, y_pred,
                        average="macro",
                        zero_division=0),
    "F1": lambda y_true, y_pred:
        f1_score(y_true, y_pred, average="macro"),
    "ACC": accuracy_score
}

results = {
    metric: {
        "train": {m: [] for m in resampling_methods},
        "test": {m: [] for m in resampling_methods}
    }
    for metric in metrics
}

# --- BUCLE PRINCIPAL ---
for run in range(n_runs):

    run_seed = random_state_base + run

    for method_name, (sampler, clf) in resampling_methods.items():

        print(f"Run {run+1}/{n_runs} - {method_name}")

        if hasattr(sampler, "random_state"):
            sampler.random_state = run_seed
        clf.random_state = run_seed

        try:

            # 1 Resampling
            X_res, y_res = sampler.fit_resample(X_train, y_train)

            # 2 Feature Selection
            selector = SelectKBest(f_classif, k=n_best_features)
            X_res_sel = selector.fit_transform(X_res, y_res)
            X_test_sel = selector.transform(X_test)

            selected_features = X_res.columns[selector.get_support()]

            # 3 Entrenamiento
            clf.fit(X_res_sel, y_res)

            train_pred = clf.predict(X_res_sel)
            test_pred = clf.predict(X_test_sel)

            y_train_eval = y_res
            y_test_eval = y_test

        except Exception as e:

            print(f"Error en {method_name} run {run+1}: {e}")
            print("-> Aplicando fallback baseline")

            clf.fit(X_train_selected_init, y_train)

            train_pred = clf.predict(X_train_selected_init)
            test_pred = clf.predict(X_test_selected_init)

            y_train_eval = y_train
            y_test_eval = y_test

            selected_features = selected_features_init

        # --- Guardado métricas (siempre) ---
        for metric_name, metric_func in metrics.items():

            results[metric_name]["train"][method_name].append(
                metric_func(y_train_eval, train_pred)
            )

            results[metric_name]["test"][method_name].append(
                metric_func(y_test_eval, test_pred)
            )

        # --- JSON metadata ---
        run_info = {
            "dataset": file_name,
            "method": method_name,
            "workflow": path_workflow,
            "run_index": run + 1,
            "random_seed": run_seed,
            "n_best_features": n_best_features,
            "selected_features": list(selected_features),
            "train_BA": results["BA"]["train"][method_name][-1],
            "test_BA": results["BA"]["test"][method_name][-1]
        }

        json_name = f"{file_name}_{method_name}-{path_workflow}_run{run+1}.json"

        with open(os.path.join(metadata_dir, json_name), "w") as f:
            json.dump(run_info, f, indent=4)

    # =============================
    # GUARDADO FINAL UNIFICADO
    # =============================
for metric_name in metrics:

    for split in ["train", "test"]:

        df_metric = pd.DataFrame(results[metric_name][split])

        df_metric.loc["max"] = df_metric.max()
        df_metric.loc["min"] = df_metric.min()
        df_metric.loc["mean"] = df_metric.mean()
        df_metric.loc["std"] = df_metric.std()

        csv_name = (
            f"{split}_{metric_name}_{file_name}"
            f"-resampling_FR-{n_best_features}.csv"
        )

        df_metric.to_csv(
            os.path.join(current_output_dir, csv_name),
            index=False
        )

print(" Ejecución 1.3 completada correctamente.")