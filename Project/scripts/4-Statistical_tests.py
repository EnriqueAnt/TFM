#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import argparse
import json
from scipy.stats import wilcoxon 
from statistics import mean
import sys

# --- ARGUMENTOS ---
parser = argparse.ArgumentParser(description='Ejecutar Tests Estadísticos (El Juez)')
parser.add_argument('--dataset', type=str, required=True, help='Nombre del dataset (ej: braaksc)')
parser.add_argument('--k', type=int, required=True, help='K features')
parser.add_argument('--input_dir', type=str, default='results_ml_classification', help='Raíz de resultados')
parser.add_argument('--output_root', type=str, default='results_ml_classification', help='Raíz de salida')

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


# --- CONFIGURACIÓN DE RUTAS DINÁMICAS ---
experiment_dir = os.path.join(args.input_dir, args.dataset, f"k_{args.k}")
stats_output_dir = os.path.join(experiment_dir, "statistics")
os.makedirs(stats_output_dir, exist_ok=True)

best_analysis_dir = os.path.join(experiment_dir, "best_analysis")
os.makedirs(best_analysis_dir, exist_ok=True)

print(f"--- INICIO JUEZ (1.4) ---")
print(f"Dataset: {args.dataset} | K: {args.k}")
print(f"Leyendo de: {experiment_dir}")
print(f"Guardando stats en: {stats_output_dir}")

# ==========================================
# 1. DEFINICIONES Y CLASES
# ==========================================

def mejor(x,y): return x > y # maximización
def peor(x,y): return x < y

class Ranking:
    def __init__(self, name):
        self.name = name
        self.wins = 0
        self.losses = 0
    def __lt__(self,x):
        return (self.wins-self.losses) < (x.wins-x.losses)
    def __str__(self):
        return f"{self.name:<35} {self.wins:>5} {self.losses:>5} {self.wins - self.losses:>5}"

def CalculateWinsLossesMatrixMean(scores):
    labels = scores.columns.values
    nScores = len(labels)
    WinLossMatriz = np.zeros((nScores,nScores)) 
    for i in range(nScores-1):            
        score_i = scores.iloc[:,i].values
        for j in range(i+1,nScores):
            score_j = scores.iloc[:,j].values
            if mejor(mean(score_i),mean(score_j)):
                WinLossMatriz[i,j] = 1
                WinLossMatriz[j,i] = -1
            if peor(mean(score_i),mean(score_j)):
                WinLossMatriz[i,j] = -1
                WinLossMatriz[j,i] = 1
    return WinLossMatriz

def CalculateWinsLossesAmount(WinsLossesMatriz,labels):
    nScores = len(labels)
    WinsLossesAmount = [Ranking(scoreName) for scoreName in labels]
    for i in range(nScores-1):            
        for j in range(i+1,nScores):
            if WinsLossesMatriz[i,j]==1:            
                WinsLossesAmount[i].wins += 1
                WinsLossesAmount[j].losses += 1
            if WinsLossesMatriz[i,j]==-1:            
                WinsLossesAmount[j].wins += 1
                WinsLossesAmount[i].losses += 1
    return WinsLossesAmount

def CalculateWinsLossesMatrixStat(scores, stat):
    labels = scores.columns.values
    nScores = len(labels)
    WinLossMatriz = np.zeros((nScores,nScores)) 
    pValues = np.zeros((nScores,nScores)) 
    for i in range(nScores-1):            
        score_i = scores.iloc[:,i].values
        for j in range(i+1,nScores):
            score_j = scores.iloc[:,j].values
            if not np.array_equal(score_i, score_j):
                try:
                    _, p_value = stat(score_i, score_j)
                except ValueError:
                    p_value = 1.0 
                
                pValues[i,j] = p_value
                pValues[j,i] = p_value                
                if p_value < 0.05:
                    if mejor(mean(score_i),mean(score_j)):
                        WinLossMatriz[i,j] = 1
                        WinLossMatriz[j,i] = -1                        
                    if peor(mean(score_i),mean(score_j)):
                        WinLossMatriz[i,j] = -1
                        WinLossMatriz[j,i] = 1                        
    return WinLossMatriz, pValues

def PrintMatriz(WinLossMatriz, labels, f_out):
    f_out.write("\nMatriz Wins-Losses:\n")
    n = len(labels)
    col_width = 25
    header = " " * col_width + "".join([f"{l[:20]:>{col_width}}" for l in labels])
    f_out.write(header + "\n")
    for i in range(n):
        row = f"{labels[i][:20]:<{col_width}}"
        for j in range(n):
            if i == j: val = "-"
            elif WinLossMatriz[j,i] == 1: val = "win"
            elif WinLossMatriz[j,i] == -1: val = "loss"
            else: val = "tie"
            row += f"{val:>{col_width}}"
        f_out.write(row + "\n")

def PrintRanking(WinLoss, f_out):
    Ranking_list = sorted(WinLoss, reverse=True)    
    f_out.write(f"\n{'Ranking':<35} {'Wins':>5} {'Losses':>5} {'Diff':>5}\n")
    for r in Ranking_list:
        f_out.write(str(r) + "\n")

# ==========================================
# 2. COMBINACIÓN DE RESULTADOS
# ==========================================

metricas = ["BA", "F1", "PS", "ACC"]
file_name_base = f"genes-{args.dataset}"

rutas = {
    "R1_Baseline": f"test_{{metrica}}_{file_name_base}_k{args.k}.csv",
    "R2_FR":       f"test_{{metrica}}_{file_name_base}-FR-{args.k}.csv",
    "R3_Resampl":  f"test_{{metrica}}_{file_name_base}-resampling_FR-{args.k}.csv"
}

combined_files = [] 

print("\n--- COMBINANDO CSVs ---")
for metrica in metricas:
    dfs = []
    found_any = False
    
    for nombre_rama, pattern in rutas.items():
        filename = pattern.format(metrica=metrica)
        file_path = os.path.join(experiment_dir, filename)
        
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            suffix = f"-{nombre_rama}"
            df.columns = [f"{col}{suffix}" for col in df.columns]
            dfs.append(df)
            found_any = True
        else:
            print(f"    [Aviso] Faltante: {filename}")

    if found_any:
        df_combinado_full = pd.concat(dfs, axis=1)
        df_combinado = df_combinado_full.iloc[:-4, :] 
        
        out_name = f"combined_{metrica}_runs.csv"
        path_out = os.path.join(stats_output_dir, out_name)
        df_combinado.to_csv(path_out, index=False)
        combined_files.append((metrica, path_out))
        print(f" Generado: {out_name}")

# ==========================================
# 3. EJECUCIÓN DE TESTS ESTADÍSTICOS
# ==========================================

report_path = os.path.join(stats_output_dir, "statistical_report.txt")

with open(report_path, 'w', encoding='utf-8') as f_log:
    f_log.write("="*60 + "\n")
    f_log.write(f"REPORTE ESTADÍSTICO - DATASET: {args.dataset} - K: {args.k}\n")
    f_log.write("="*60 + "\n")

    for metrica, file_path in combined_files:
        f_log.write(f"\n\n>>> ANALIZANDO MÉTRICA: {metrica} <<<\n")
        f_log.write("-" * 40 + "\n")
        
        scores = pd.read_csv(file_path)
        
        # 1. Ranking Medias
        f_log.write("\n--- Ranking de Medias (Wins-Losses) ---\n")
        WinLossMatrizMean = CalculateWinsLossesMatrixMean(scores)
        WinsLossesAmountMean = CalculateWinsLossesAmount(WinLossMatrizMean, scores.columns.values)
        PrintMatriz(WinLossMatrizMean, scores.columns.values, f_log)
        PrintRanking(WinsLossesAmountMean, f_log)
        
        # 2. Ranking Wilcoxon
        f_log.write("\n--- Ranking Wilcoxon (p < 0.05) ---\n")
        WinLossMatrizStat, pValues = CalculateWinsLossesMatrixStat(scores, wilcoxon)
        WinsLossesAmountStat = CalculateWinsLossesAmount(WinLossMatrizStat, scores.columns.values)
        PrintMatriz(WinLossMatrizStat, scores.columns.values, f_log)
        PrintRanking(WinsLossesAmountStat, f_log)

print(f" Reporte estadístico guardado en: {report_path}")

# ==========================================
# 4. DETERMINACIÓN DEL GANADOR
# ==========================================

print("\n--- DETERMINANDO GANADOR GLOBAL ---")

resultados_signed = {}
resultados_mean = {}
metricas_resumen = ["BA", "F1", "PS"] 

for metrica in metricas_resumen:
    path = next((p for m, p in combined_files if m == metrica), None)
    
    if path:
        scores = pd.read_csv(path)
        labels = scores.columns.values
        
        matrix_sr, _ = CalculateWinsLossesMatrixStat(scores, wilcoxon)
        signed_summary = CalculateWinsLossesAmount(matrix_sr, labels)
        resultados_signed[metrica] = {r.name: r.wins - r.losses for r in signed_summary}
        
        matrix_mean = CalculateWinsLossesMatrixMean(scores)
        mean_summary = CalculateWinsLossesAmount(matrix_mean, labels)
        resultados_mean[metrica] = {r.name: r.wins - r.losses for r in mean_summary}

df_signed = pd.DataFrame(resultados_signed).fillna(0)
df_mean = pd.DataFrame(resultados_mean).fillna(0)

# --- CÁLCULO DE TOTALES ---
df_signed["Total"] = df_signed.sum(axis=1)
df_mean["Total"] = df_mean.sum(axis=1)

# Agregar columnas temporales para el desempate final
df_signed["Total_Mean"] = df_mean["Total"]
df_mean["Total_SR"] = df_signed["Total"]

# Asegurar que BA y F1 existan en caso de que alguna ejecución parcial no las haya calculado
sort_cols_sr = ["Total"]
sort_cols_mean = ["Total"]

if "BA" in df_signed.columns:
    sort_cols_sr.append("BA")
    sort_cols_mean.append("BA")
if "F1" in df_signed.columns:
    sort_cols_sr.append("F1")
    sort_cols_mean.append("F1")

sort_cols_sr.append("Total_Mean")
sort_cols_mean.append("Total_SR")

# --- ORDENAR Y GUARDAR SUMMARIES CSV (Visualización de mayor a menor) ---
df_signed_sorted = df_signed.sort_values(by=sort_cols_sr, ascending=[False]*len(sort_cols_sr))
df_mean_sorted = df_mean.sort_values(by=sort_cols_mean, ascending=[False]*len(sort_cols_mean))

# Eliminar las columnas cruzadas antes de guardar
df_signed_sorted = df_signed_sorted.drop(columns=["Total_Mean"])
df_mean_sorted = df_mean_sorted.drop(columns=["Total_SR"])
df_signed = df_signed.drop(columns=["Total_Mean"])
df_mean = df_mean.drop(columns=["Total_SR"])

path_summary_sr = os.path.join(best_analysis_dir, "summary_signed_rank.csv")
path_summary_mean = os.path.join(best_analysis_dir, "summary_mean.csv")

df_signed_sorted.to_csv(path_summary_sr)
df_mean_sorted.to_csv(path_summary_mean)

print(f" Resúmenes CSV generados y ordenados en best_analysis.")

winner_final = None
best_run_num = -1
criterio = ""

if not df_signed.empty:
    # Mantener el orden interno para la selección lógica del ganador
    df_signed = df_signed_sorted
    df_mean = df_mean_sorted

    winner_sr = df_signed.iloc[0]
    winner_mean = df_mean.iloc[0]
    
    diff_mean_rel = df_mean.loc[winner_mean.name, "Total"] - df_mean.loc[winner_sr.name, "Total"]
    diff_sr_rel = winner_sr["Total"] - df_signed.loc[winner_mean.name, "Total"]

    criterio = "Signed-Rank"
    winner_final = winner_sr

    if winner_mean.name != winner_sr.name:
        if diff_mean_rel > 3 * diff_sr_rel:
            winner_final = winner_mean
            criterio = "Mean (Diff-Mean > 3x Diff-SR)"

    print(f" Ganador Final: {winner_final.name}")
    print(f"   Criterio: {criterio}")

    # ==========================================
    # 5. RECONSTRUCCIÓN DEL MEJOR RUN Y JSON
    # ==========================================
    
    metrica_ref = "BA"
    path_ref = next((p for m, p in combined_files if m == metrica_ref), None)
    
    if path_ref:
        scores_ref = pd.read_csv(path_ref)
        if winner_final.name in scores_ref.columns:
            run_values = scores_ref[winner_final.name].values
            best_run_index = int(np.argmax(run_values)) # 0-9
            best_run_num = best_run_index + 1 # 1-10
            
            print(f" Mejor Run: #{best_run_num} (Score: {run_values[best_run_index]:.4f})")
            
            # Buscar JSON (Respetando la carpeta runs_metadata local)
            json_runs_dir = os.path.join(experiment_dir, "runs_metadata")
            target_method_part = winner_final.name.split('-')[0]
            target_branch_part = winner_final.name.split('-')[1].split('_')[0] # Extraer R1, R2 o R3
            
            found_json_data = None
            
            # CASO R2 o R3: Buscar en JSONs existentes
            if "R1" not in winner_final.name and os.path.exists(json_runs_dir):
                possible_jsons = [f for f in os.listdir(json_runs_dir) if f.endswith(f"_run{best_run_num}.json")]
                for pj in possible_jsons:
                    if target_method_part in pj:
                        with open(os.path.join(json_runs_dir, pj), 'r') as f:
                            found_json_data = json.load(f)
                            # Forzar clave path_workflow coherente con el ganador real
                            found_json_data["path_workflow"] = target_branch_part
                            if "workflow" in found_json_data:
                                del found_json_data["workflow"]
                        print(f" JSON encontrado: {pj}")
                        break
            
            # CASO R1: Generar JSON manual, todos los genes.
            elif "R1" in winner_final.name:
                print(f" R1 Ganador: Generando metadatos de genes totales...")
                dataset_original_path = os.path.join("data", f"genes-{args.dataset}.csv")
                all_genes = []
                
                if os.path.exists(dataset_original_path):
                    df_all = pd.read_csv(dataset_original_path, nrows=0)
                    all_genes = [c for c in df_all.columns if c not in ['target', 'sample_id', 'Class', 'label', 'Unnamed: 0']]
                
                seed_f1 = 42 + best_run_index
                
                found_json_data = {
                    "dataset": f"genes-{args.dataset}",
                    "method": target_method_part,
                    "path_workflow": "R1",
                    "run_index": best_run_index,
                    "random_seed": seed_f1,
                    "n_best": 0,
                    "n_best_features": args.k,
                    "n_best_features_full": len(all_genes),
                    "selected_features": all_genes
                }

            # Guardar el archivo final si tenemos datos
            if found_json_data:
                with open(os.path.join(best_analysis_dir, "best_model_info.txt"), "w") as f:
                    json.dump(found_json_data, f, indent=4)

            
    # =======================================================
    # 6. GENERACIÓN DEL ARCHIVO .PY (BASE DE DATOS GANADOR)
    # =======================================================
    print("\n--- GENERANDO FICHA TÉCNICA (.PY) ---")
    
    # Extraer métricas específicas del mejor run para todas las métricas disponibles
    metrics_run_values = {} # diccionario {BA: {train: 0.9, test: 0.8}, ...}
    
    for m in ["BA", "F1", "PS", "ACC"]:
        path_test = next((p for met, p in combined_files if met == m), None)
        
        parts = winner_final.name.split('-')
        method_raw = parts[0] # SMOTE_RF
        branch_raw = parts[1].split('_')[0] # R1, R2 o R3
        
        val_test = 0.0
        if path_test:
            sdf = pd.read_csv(path_test)
            if winner_final.name in sdf.columns:
                val_test = sdf.iloc[best_run_index][winner_final.name]
        
        file_pattern = ""
        if "R1" in winner_final.name: file_pattern = f"train_{{m}}_{file_name_base}_k{args.k}.csv"
        elif "R2" in winner_final.name: file_pattern = f"train_{{m}}_{file_name_base}-FR-{args.k}.csv"
        elif "R3" in winner_final.name: file_pattern = f"train_{{m}}_{file_name_base}-resampling_FR-{args.k}.csv"
        
        val_train = 0.0
        try:
            train_file = os.path.join(experiment_dir, file_pattern.format(m=m))
            if os.path.exists(train_file):
                df_train = pd.read_csv(train_file)
                if method_raw in df_train.columns:
                    val_train = df_train.iloc[best_run_index][method_raw]
        except Exception as e:
            pass 
            
        metrics_run_values[m] = {"train": val_train, "test": val_test}

    py_filename = f"{args.dataset}_best_dataset_k{args.k}.py"
    py_path = os.path.join(best_analysis_dir, py_filename)
    
    path_data_ref = experiment_dir.replace("\\", "/")
    
    py_content = f"""# Información del Mejor Dataset y Modelo Automático
# Generado por script 4-Statistical_tests.py

dataset_name = '{args.dataset}'
k_value = {args.k}

# Rutas de referencia (Carpeta del experimento)
experiment_path = '{path_data_ref}'
# Nota: Para R3, los datasets cambian por run. Para R1/R2 son fijos.
# Ruta aproximada de los datos procesados:
train_dataset_path_template = '{path_data_ref}/{file_name_base}_train.csv' 

winner_methodology = '{winner_final.name}'
winner_method = '{method_raw}'
winner_branch = '{branch_raw}'

best_run_index = '{best_run_index}'
Run_info = 'Run #{best_run_num} (equivale a {best_run_num}/10)'

# Métricas del ganador (Run #{best_run_num})
winner_train_ba = {metrics_run_values['BA']['train']}
winner_test_ba  = {metrics_run_values['BA']['test']}

winner_train_ps = {metrics_run_values['PS']['train']}
winner_test_ps  = {metrics_run_values['PS']['test']}

winner_train_f1 = {metrics_run_values['F1']['train']}
winner_test_f1  = {metrics_run_values['F1']['test']}

winner_train_acc = {metrics_run_values['ACC']['train']}
winner_test_acc  = {metrics_run_values['ACC']['test']}
"""

    with open(py_path, "w", encoding="utf-8") as f:
        f.write(py_content)
    
    print(f" Archivo Python generado: {py_path}")

else:
    print(" No se encontraron resultados suficientes para calcular ganador.")

print("Ejecución del Juez finalizada.")