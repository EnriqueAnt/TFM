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

import random
import numpy as np
import pandas as pd
import argparse

# 2. FUNCIÓN DE REPRODUCIBILIDAD GLOBAL
def set_reproducibility(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)

# Fijamos la semilla base coincidente con el resto del proyecto
set_reproducibility(42)

# 3. IMPORTACIONES DE SKLEARN E IMBLEARN
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import balanced_accuracy_score, precision_score, f1_score, accuracy_score
from imblearn.ensemble import (
    BalancedRandomForestClassifier, 
    EasyEnsembleClassifier, 
    RUSBoostClassifier, 
    BalancedBaggingClassifier
)
from sklearn.tree import DecisionTreeClassifier

# --- BLOQUE DE ARGUMENTOS ---
parser = argparse.ArgumentParser(description='Ejecutar experimento 1.1 COMPARATIVA DE CLASIFICADORES (SIN SELECCIÓN DE K)')
parser.add_argument('--dataset', type=str, required=True, help='Nombre del dataset')
parser.add_argument('--k', type=int, required=True, help='Se recibe por compatibilidad pero SE IGNORA en 1.1')
parser.add_argument('--input_dir', type=str, default='.', help='Directorio datos entrada')
parser.add_argument('--output_root', type=str, default='results_ml_classification', help='Directorio resultados')

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


# Variables
file_name = f"genes-{args.dataset}"
# NOTA: En 1.1 ignoramos K para la selección, 
# pero usamos el valor de K para nombrar la carpeta de salida correctamente.
k_label = args.k 

# Estructura de carpetas
current_output_dir = os.path.join(args.output_root, args.dataset, f"k_{k_label}")
os.makedirs(current_output_dir, exist_ok=True)

print(f"--- INICIANDO SCRIPT 1.1 COMPARATIVA DE CLASIFICADORES (SIN SELECCIÓN DE K)---")
print(f"Dataset: {args.dataset}")
print(f"Nota: Se usarán TODAS las características (K={k_label} recibido pero ignorado en filtrado)")
print(f"Guardando en: {current_output_dir}")

# Parámetros generales
random_state_base = 42
n_runs = 10

# Carga de datos
input_path = os.path.join(args.input_dir, f"{file_name}.csv")
print(f"Leyendo archivo: {input_path}")
df = pd.read_csv(input_path)

# Preparar X e y
X = df.iloc[:, :-1]
y = df.iloc[:, -1]
X = pd.get_dummies(X)

# Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=random_state_base)

# --- SIN SELECCIÓN DE CARACTERÍSTICAS (CLASIFICADORES) ---
print(f"Usando todas las características disponibles: {X_train.shape[1]}")
X_train_selected = X_train.copy()
X_test_selected = X_test.copy()

print(f"Dimensiones Train: {X_train_selected.shape}")

# Guardar datasets procesados (Raw)
train_selected_df = pd.DataFrame(X_train_selected)
train_selected_df['target'] = y_train.values

test_selected_df = pd.DataFrame(X_test_selected)
test_selected_df['target'] = y_test.values

# Guardamos los CSVs. 
# run_all.py espera encontrar archivos que terminan en _train.csv o _test.csv
train_output_path = os.path.join(current_output_dir, f"{file_name}_train.csv")
test_output_path = os.path.join(current_output_dir, f"{file_name}_test.csv")

train_selected_df.to_csv(train_output_path, index=False)
test_selected_df.to_csv(test_output_path, index=False)

# Definición de clasificadores
classifiers = {
    "RF_None": RandomForestClassifier(class_weight=None), 
    "RF_Balanced": RandomForestClassifier(class_weight="balanced"),
    "RF_BalancedSub": RandomForestClassifier(class_weight="balanced_subsample"),
    "BRF_None": BalancedRandomForestClassifier(class_weight=None, sampling_strategy="all", replacement=True, bootstrap=False), 
    "BRF_Balanced": BalancedRandomForestClassifier(class_weight="balanced", sampling_strategy="all", replacement=True, bootstrap=False), 
    "BRF_BalancedSub": BalancedRandomForestClassifier(class_weight="balanced_subsample", sampling_strategy="all", replacement=True, bootstrap=False),
    "SVC_Balanced": SVC(probability=True, class_weight="balanced"),    
    "HGB_Balanced": HistGradientBoostingClassifier(class_weight="balanced"),
    "EEC": EasyEnsembleClassifier(),
    "RBC": RUSBoostClassifier(
        estimator=DecisionTreeClassifier(max_depth=3), 
        n_estimators=50,                               
        learning_rate=0.1,                             
        sampling_strategy='auto',
        random_state=42                                
        ),
    "BBC": BalancedBaggingClassifier(),
}
# Aplicar regularización a los clasificadores
def apply_regularization(clf):
    
    # Random Forest (RF_None, RF_Balanced, RF_BalancedSub) (sin sampling externo explícito)
    if isinstance(clf, RandomForestClassifier): 
        clf.set_params(
            n_estimators=100,
            max_depth=10,
            min_samples_leaf=5,
            min_samples_split=10,
            max_features='sqrt'
        )

    # Balanced Random Forest (BRF_None, BRF_Balanced, BRF_BalancedSub) (con re-sampling interno)
    elif isinstance(clf, BalancedRandomForestClassifier):
        clf.set_params(
            n_estimators=100,
            max_depth=10,
            min_samples_leaf=5,
            min_samples_split=10
        )

    # SVM (SVC_Balanced)
    elif isinstance(clf, SVC):
        clf.set_params(
            C=1.0,
            kernel='linear'
        )

    # HistGradientBoosting (HGB_Balanced)
    elif isinstance(clf, HistGradientBoostingClassifier):
        clf.set_params(
            max_depth=6,
            learning_rate=0.05,
            max_iter=100
        )

    # EasyEnsemble (EEC)
    elif isinstance(clf, EasyEnsembleClassifier):
        clf.set_params(
            n_estimators=10
        )

    # RUSBoost (RBC)
    elif isinstance(clf, RUSBoostClassifier):
        clf.set_params(
            estimator=DecisionTreeClassifier(max_depth=3)
        )

    # BalancedBagging (BBC)
    elif isinstance(clf, BalancedBaggingClassifier):
        clf.set_params(
            estimator=DecisionTreeClassifier(max_depth=5),
            n_estimators=10
        )

    return clf

# Inicializar almacenamiento
metrics = ["ba", "ps", "f1", "acc"]
train_scores = {m: {clf: [] for clf in classifiers} for m in metrics}
test_scores = {m: {clf: [] for clf in classifiers} for m in metrics}

# Bucle de entrenamiento
for run in range(n_runs):
    run_seed = random_state_base + run
    
    for clf_name, clf in classifiers.items():
        print(f"Entrenando {clf_name} - Run {run+1}")
        clf = apply_regularization(clf)
        clf.random_state = run_seed
        
        try:
            clf.fit(X_train_selected, y_train)
            
            train_pred = clf.predict(X_train_selected)
            test_pred = clf.predict(X_test_selected)
            
            # Métricas
            train_scores["ba"][clf_name].append(balanced_accuracy_score(y_train, train_pred))
            test_scores["ba"][clf_name].append(balanced_accuracy_score(y_test, test_pred))
            
            train_scores["ps"][clf_name].append(precision_score(y_train, train_pred, average='macro', zero_division=np.nan))
            test_scores["ps"][clf_name].append(precision_score(y_test, test_pred, average='macro', zero_division=np.nan))
            
            train_scores["f1"][clf_name].append(f1_score(y_train, train_pred, average='macro'))
            test_scores["f1"][clf_name].append(f1_score(y_test, test_pred, average='macro'))
            
            train_scores["acc"][clf_name].append(accuracy_score(y_train, train_pred))
            test_scores["acc"][clf_name].append(accuracy_score(y_test, test_pred))
            
        except Exception as e:
            print(f" AVISO: Fallo en {clf_name} Run {run+1}. Error: {e}")
            # Penalización
            for m in metrics:
                train_scores[m][clf_name].append(0.0)
                test_scores[m][clf_name].append(0.0)

# Guardar Resultados
metric_names_map = {"ba": "BA", "ps": "PS", "f1": "F1", "acc": "ACC"}

for m_key, m_name in metric_names_map.items():
    # Train
    df_train = pd.DataFrame(train_scores[m_key])
    df_train.loc["max"] = df_train.max()
    df_train.loc["min"] = df_train.min()
    df_train.loc["mean"] = df_train.mean()
    df_train.loc["std"] = df_train.std()
    
    path_train = os.path.join(current_output_dir, f"train_{m_name}_{file_name}_k{k_label}.csv")
    df_train.to_csv(path_train, index=False)
    
    # Test
    df_test = pd.DataFrame(test_scores[m_key])
    df_test.loc["max"] = df_test.max()
    df_test.loc["min"] = df_test.min()
    df_test.loc["mean"] = df_test.mean()
    df_test.loc["std"] = df_test.std()
    
    path_test = os.path.join(current_output_dir, f"test_{m_name}_{file_name}_k{k_label}.csv")
    df_test.to_csv(path_test, index=False)

print("Ejecución 1.1 Baseline & Specialized Classifier completada correctamente.")