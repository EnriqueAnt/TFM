#!/usr/bin/env python
# coding: utf-8

import os

# 1. CONFIGURACIÓN DE ENTORNO (Crítico: Antes de importar numpy o sklearn)
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
from collections import Counter

# 2. FUNCIÓN DE REPRODUCIBILIDAD GLOBAL
def set_reproducibility(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)

# Fijamos la semilla base
set_reproducibility(42)

# 3. IMPORTACIONES RESTANTES
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, precision_score, f1_score, accuracy_score
from sklearn.cluster import MiniBatchKMeans, KMeans

from imblearn.under_sampling import (
    ClusterCentroids, CondensedNearestNeighbour, EditedNearestNeighbours,
    RepeatedEditedNearestNeighbours, AllKNN, NearMiss, NeighbourhoodCleaningRule, 
    OneSidedSelection, RandomUnderSampler, TomekLinks
)
from imblearn.over_sampling import (
    RandomOverSampler, SMOTE, BorderlineSMOTE, KMeansSMOTE
)
from imblearn.combine import SMOTEENN, SMOTETomek

# --- BLOQUE DE ARGUMENTOS ---
parser = argparse.ArgumentParser(description='Ejecutar experimento 1.2 Resampling con RF')
parser.add_argument('--dataset', type=str, required=True, help='Nombre del dataset (ej: braaksc)')
parser.add_argument('--k', type=int, required=True, help='Número de características K')
parser.add_argument('--input_dir', type=str, default='.', help='Directorio datos entrada')
parser.add_argument('--output_root', type=str, default='results_ml_classification', help='Directorio raiz resultados')

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
path_workflow = "R2"
random_state_base = 42
n_runs = 10

# Carpetas
current_output_dir = os.path.join(args.output_root, args.dataset, f"k_{n_best_features}")
os.makedirs(current_output_dir, exist_ok=True)
metadata_dir = os.path.join(current_output_dir, "runs_metadata")
os.makedirs(metadata_dir, exist_ok=True)

print(f"--- INICIO SCRIPT 1.2 ---")
print(f"Dataset: {file_name} | K: {n_best_features}")

# --- CARGA DE DATOS ---
input_path = os.path.join(args.input_dir, f"{file_name}.csv")
df = pd.read_csv(input_path)

# Split
X = df.iloc[:, :-1]
y = df.iloc[:, -1]
X = pd.get_dummies(X)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=random_state_base
)

# --- SELECCIÓN DE CARACTERÍSTICAS ---
print(f"Seleccionando {n_best_features} mejores características...")
selector = SelectKBest(score_func=f_classif, k=n_best_features)
X_train_selected = selector.fit_transform(X_train, y_train)
X_test_selected = selector.transform(X_test)

selected_features = X_train.columns[selector.get_support()]

# Guardar datasets
train_selected_df = pd.DataFrame(X_train_selected, columns=selected_features)
train_selected_df['target'] = y_train.values
test_selected_df = pd.DataFrame(X_test_selected, columns=selected_features)
test_selected_df['target'] = y_test.values

train_selected_df.to_csv(os.path.join(current_output_dir, f"{file_name}-FR-{n_best_features}_train.csv"), index=False)
test_selected_df.to_csv(os.path.join(current_output_dir, f"{file_name}-FR-{n_best_features}_test.csv"), index=False)

# --- LÓGICA CONDICIONAL DE PARÁMETROS (COGDX vs RESTO) ---
if args.dataset == 'cogdx':
    print("Configuración ESPECÍFICA para COGDX aplicada (k=2, clusters=1)")
    k_neighbors_val = 2
    m_neighbors_val = 2 
    n_clusters_val = 1
    balance_threshold_val = 0.001
else:
    print(f"Configuración ESTÁNDAR aplicada para {args.dataset} (k=3, clusters=2)")
    k_neighbors_val = 3
    m_neighbors_val = 3
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
    
    "BorderlineSMOTE": BorderlineSMOTE(
        k_neighbors=k_neighbors_val, 
        m_neighbors=m_neighbors_val
    ),
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

# --- RF INDEPENDIENTE POR MÉTODO ---
resampling_methods = {
    f"{name}_RF": (
        method,
        RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            min_samples_leaf=5,
            min_samples_split=10,
            max_features='sqrt',
            class_weight=None
        )
    )
    for name, method in base_methods.items()
}

# Inicializar almacenamiento
metrics = ["ba", "ps", "f1", "acc"]
results = {m: {"train": {method: [] for method in resampling_methods}, 
               "test": {method: [] for method in resampling_methods}} for m in metrics}

# --- BUCLE DE ENTRENAMIENTO ---
for run in range(n_runs):
    run_seed = random_state_base + run

    for method_name, (sampler, clf) in resampling_methods.items():
        if hasattr(sampler, "random_state"):
            sampler.random_state = run_seed        
        clf.random_state = run_seed

        try:
            # INTENTO NORMAL: Resampling + Fit
            X_resampled, y_resampled = sampler.fit_resample(X_train_selected, y_train)
            clf.fit(X_resampled, y_resampled)
            
            train_pred = clf.predict(X_resampled)
            test_pred = clf.predict(X_test_selected)

            # Métricas
            results["ba"]["train"][method_name].append(balanced_accuracy_score(y_resampled, train_pred))
            results["ba"]["test"][method_name].append(balanced_accuracy_score(y_test, test_pred))
            results["ps"]["train"][method_name].append(precision_score(y_resampled, train_pred, average='macro', zero_division=np.nan))
            results["ps"]["test"][method_name].append(precision_score(y_test, test_pred, average='macro', zero_division=np.nan))
            results["f1"]["train"][method_name].append(f1_score(y_resampled, train_pred, average='macro'))
            results["f1"]["test"][method_name].append(f1_score(y_test, test_pred, average='macro'))
            results["acc"]["train"][method_name].append(accuracy_score(y_resampled, train_pred))
            results["acc"]["test"][method_name].append(accuracy_score(y_test, test_pred))

            # JSON Metadata
            run_info = {
                "dataset": file_name,
                "method": method_name,
                "path_workflow": path_workflow,
                "run_index": run + 1,
                "random_seed": run_seed,
                "n_best_features": n_best_features,
                "selected_features": list(selected_features),
                "train_balanced_accuracy": results["ba"]["train"][method_name][-1],
                "test_balanced_accuracy": results["ba"]["test"][method_name][-1]
            }

        except Exception as e:
            print(f" Error en {method_name} (Run {run+1}): {e} -> Usando FALLBACK (Sin Resampling)")
            
            # FALLBACK: Entrenar con datos originales seleccionados (sin balancear)
            # Esto evita que las listas queden con un elemento menos y rompan el Excel final
            clf.fit(X_train_selected, y_train)
            
            train_pred = clf.predict(X_train_selected)
            test_pred = clf.predict(X_test_selected)
            
            results["ba"]["train"][method_name].append(balanced_accuracy_score(y_train, train_pred))
            results["ba"]["test"][method_name].append(balanced_accuracy_score(y_test, test_pred))
            results["ps"]["train"][method_name].append(precision_score(y_train, train_pred, average='macro', zero_division=np.nan))
            results["ps"]["test"][method_name].append(precision_score(y_test, test_pred, average='macro', zero_division=np.nan))
            results["f1"]["train"][method_name].append(f1_score(y_train, train_pred, average='macro'))
            results["f1"]["test"][method_name].append(f1_score(y_test, test_pred, average='macro'))
            results["acc"]["train"][method_name].append(accuracy_score(y_train, train_pred))
            results["acc"]["test"][method_name].append(accuracy_score(y_test, test_pred))

            run_info = {
                "dataset": file_name,
                "method": method_name,
                "path_workflow": path_workflow,
                "run_index": run + 1,
                "error_log": str(e),
                "fallback_applied": True,
                "train_balanced_accuracy": results["ba"]["train"][method_name][-1],
                "test_balanced_accuracy": results["ba"]["test"][method_name][-1]
            }

        # Guardar JSON (sea normal o fallback)
        json_filename = f"{file_name}_{method_name}-{path_workflow}_run{run+1}.json"
        with open(os.path.join(metadata_dir, json_filename), "w", encoding="utf-8") as f:
            json.dump(run_info, f, indent=4)

# Guardar Resultados Finales
metric_map = {"ba": "BA", "ps": "PS", "f1": "F1", "acc": "ACC"}

for m_key, m_name in metric_map.items():
    for mode in ["train", "test"]:
        df_results = pd.DataFrame(results[m_key][mode])
        df_results.loc["max"] = df_results.max()
        df_results.loc["min"] = df_results.min()
        df_results.loc["mean"] = df_results.mean()
        df_results.loc["std"] = df_results.std()

        csv_filename = f"{mode}_{m_name}_{file_name}-FR-{n_best_features}.csv"
        df_results.to_csv(os.path.join(current_output_dir, csv_filename), index=False)

print("Ejecución 1.2 completada correctamente.")