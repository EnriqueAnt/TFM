import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Configuración
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results_ml_classification")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "parsimony_results")

K_VALUES = [25, 50, 75, 100, 500, 1000, 1500, 2500]
DATASETS = ["braaksc", "ceradsc", "cogdx"]
PARSIMONY_THRESHOLD = 0.01

def get_performance_metrics(dataset, k):
    """Extrae las métricas de rendimiento y calcula la puntuación ponderada."""
    py_path = os.path.join(RESULTS_DIR, dataset, f"k_{k}", "best_analysis", f"{dataset}_best_dataset_k{k}.py")
    if not os.path.exists(py_path):
        return None
    try:
        with open(py_path, 'r', encoding='utf-8') as f:
            content = f.read()
            ba = float(re.search(r"winner_test_ba\s*=\s*([\d\.]+)", content).group(1))
            f1 = float(re.search(r"winner_test_f1\s*=\s*([\d\.]+)", content).group(1))
            ps = float(re.search(r"winner_test_ps\s*=\s*([\d\.]+)", content).group(1))
        
        # Puntuación ponderada: 50% BA, 35% F1, 15% Precisión
        score = (0.50 * ba + 0.35 * f1 + 0.15 * ps)
        return {"BA": ba, "F1": f1, "PS": ps, "Score": score}
    except Exception:
        return None

def find_knee_point(ks, scores):
    """
    Encuentra el 'Knee Point' (codo) usando la distancia perpendicular positiva máxima 
    desde la línea de tendencia en el espacio logarítmico.
    """
    n_points = len(ks)
    log_ks = np.log10(ks)
    coords = np.column_stack((log_ks, scores))
    
    first_point = coords[0]
    last_point = coords[-1]
    
    line_vec = last_point - first_point
    line_vec_norm = line_vec / np.linalg.norm(line_vec)
    
    # Vector normal (rotación de 90 grados) para medir la distancia con signo
    normal_vec = np.array([-line_vec_norm[1], line_vec_norm[0]])
    
    vec_from_first = coords - first_point
    signed_distances = np.dot(vec_from_first, normal_vec)
    
    # Retorna K con la mayor "curva" positiva
    return ks[np.argmax(signed_distances)]

def run_parsimony_validation():
    """Función principal de ejecución para la auditoría de parsimonia."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    summary_data = []
    
    print(f"{'Dataset':<10} | {'Max Score':<10} | {'Best K':<8} | {'Ockham K':<10} | {'Knee K':<8} | {'Status'}")
    print("-" * 70)

    for ds in DATASETS:
        dataset_results = []
        for k in K_VALUES:
            metrics = get_performance_metrics(ds, k)
            if metrics:
                metrics["K"] = k
                dataset_results.append(metrics)
        
        if not dataset_results:
            continue
            
        df = pd.DataFrame(dataset_results).sort_values("K")
        ks = df["K"].values
        scores = df["Score"].values
        
        # 1. Navaja de Ockham (Umbral del 1%)
        best_idx = np.argmax(scores)
        best_score = scores[best_idx]
        best_k = ks[best_idx]
        threshold = best_score - PARSIMONY_THRESHOLD
        
        # El K más pequeño dentro del 1% de la mejor puntuación
        ockham_k = ks[np.where(scores >= threshold)[0][0]]
        
        # 2. Punto Knee Geométrico (Codo)
        knee_k = find_knee_point(ks, scores)
        
        status = "MATCH" if ockham_k == knee_k else "DIVERGENT"
        print(f"{ds:<10} | {best_score:<10.4f} | {best_k:<8} | {ockham_k:<10} | {knee_k:<8} | {status}")
        
        summary_data.append({
            "Dataset": ds,
            "Max_Score": best_score,
            "Raw_Best_K": best_k,
            "Ockham_K": ockham_k,
            "Knee_Point_K": knee_k,
            "Status": status
        })
        
        # --- Visualization ---
        plt.figure(figsize=(10, 6))
        plt.plot(ks, scores, marker='o', linestyle='-', color='#7f8c8d', alpha=0.6, linewidth=2, label='Weighted Score (Test)')
        
        # Reference lines
        plt.axhline(y=best_score, color='#e74c3c', linestyle='--', alpha=0.4, label=f'Max ({best_score:.4f})')
        plt.axhline(y=threshold, color='#27ae60', linestyle=':', linewidth=2, label='1% Parsimony Threshold')
        plt.fill_between(ks, threshold, best_score, color='#2ecc71', alpha=0.1)
        
        # Highlight key points
        plt.scatter([ockham_k], [scores[ks==ockham_k]], color='#27ae60', s=180, zorder=5, edgecolors='black', label=f'Ockham K ({ockham_k})')
        plt.scatter([knee_k], [scores[ks==knee_k]], color='#2980b9', s=120, marker='D', zorder=5, edgecolors='black', label=f'Knee Point K ({knee_k})')
        
        plt.title(f"Parsimony Validation Analysis: {ds.upper()}", fontsize=14, fontweight='bold')
        plt.xlabel("Number of Selected Features (K)", fontsize=12)
        plt.ylabel("Performance Score", fontsize=12)
        plt.xscale('log')
        plt.xticks(K_VALUES, K_VALUES)
        plt.legend(loc='lower right', frameon=True, shadow=True)
        plt.grid(True, which="both", ls="-", alpha=0.1)
        
        plt.savefig(os.path.join(OUTPUT_DIR, f"parsimony_plot_{ds}.png"), dpi=300, bbox_inches='tight')
        plt.close()

    # Save summary report
    pd.DataFrame(summary_data).to_csv(os.path.join(OUTPUT_DIR, "parsimony_audit_summary.csv"), index=False)
    print(f"\nAudit completed. Results stored in: {OUTPUT_DIR}")

if __name__ == "__main__":
    run_parsimony_validation()
