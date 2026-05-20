import os
import pandas as pd
import re
import shutil
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np

# ==================== CONFIGURACIÓN ====================
# Obtenemos la ruta absoluta de la carpeta donde está este script (Project/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

datasets       = ["braaksc", "ceradsc", "cogdx"]
k_values       = [25, 50, 75, 100, 500, 1000, 1500, 2500]

results_root   = os.path.join(SCRIPT_DIR, "results_ml_classification")
final_output_root = os.path.join(SCRIPT_DIR, "best_k_selection")
PARSIMONY_THRESHOLD = 0.01   # 1% tolerancia

# BA pondera más porque el dataset puede tener clases desbalanceadas
# F1 captura el equilibrio precisión/recall
# PS (Precision) es secundaria en problemas clínicos donde los falsos negativos son costosos
METRIC_WEIGHTS = {"BA": 0.50, "F1": 0.35, "PS": 0.15}

# Subcarpeta nueva para todas las visualizaciones y tablas de este script
ANALYSIS_SUBDIR = os.path.join(final_output_root, "parsimony_analysis")

# ==================== GESTIÓN SEGURA DE DIRECTORIOS ====================

def setup_directories():
    """Crea la estructura de carpetas de forma segura, sin borrar lo existente."""

    if os.path.exists(final_output_root):
        print(f"  La carpeta '{final_output_root}' ya existe. Se preservarán los datos anteriores.")
        print(f"    Solo se sobreescribirá la subcarpeta '{ANALYSIS_SUBDIR}'.\n")
    else:
        os.makedirs(final_output_root)

    # La subcarpeta de análisis sí se regenera limpia en cada ejecución
    if os.path.exists(ANALYSIS_SUBDIR):
        shutil.rmtree(ANALYSIS_SUBDIR)
    os.makedirs(ANALYSIS_SUBDIR)
    os.makedirs(os.path.join(ANALYSIS_SUBDIR, "plots"))

# ==================== FUNCIONES DE PARSEO ====================

def parse_best_dataset_py(filepath):
    """Lee el archivo .py del Juez y extrae las variables clave."""
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        match_winner = re.search(r"winner_methodology\s*=\s*['\"](.+?)['\"]", content)
        match_ba      = re.search(r"winner_test_ba\s*=\s*([\d\.]+)", content)
        match_f1      = re.search(r"winner_test_f1\s*=\s*([\d\.]+)", content)
        match_ps      = re.search(r"winner_test_ps\s*=\s*([\d\.]+)", content)

        if not all([match_winner, match_ba, match_f1, match_ps]):
            return None

        return {
            'Ganador': match_winner.group(1),
            'BA_Test': float(match_ba.group(1)),
            'F1_Test': float(match_f1.group(1)),
            'PS_Test': float(match_ps.group(1))
        }
    except Exception as e:
        print(f"      Error al parsear {filepath}: {e}")
        return None


def get_dataset_filename_pattern(dataset, k, winner_name):
    """Mapea el nombre del ganador al archivo CSV original."""
    base = f"genes-{dataset}"
    if "R1" in winner_name:
        return f"{base}_train.csv", f"{base}_test.csv"
    elif "R2" in winner_name:
        return f"{base}-FR-{k}_train.csv", f"{base}-FR-{k}_test.csv"
    elif "R3" in winner_name:
        return f"{base}-resampling_FR-{k}_train.csv", f"{base}-resampling_FR-{k}_test.csv"
    else:
        return f"{base}_train.csv", f"{base}_test.csv"

# ==================== CÁLCULO DE SCORE PONDERADO ====================

def weighted_score(ba, f1, ps):
    """Score ponderado: BA tiene más peso por posible desbalance de clases."""
    return (METRIC_WEIGHTS["BA"] * ba +
            METRIC_WEIGHTS["F1"] * f1 +
            METRIC_WEIGHTS["PS"] * ps)

# ==================== VISUALIZACIÓN: GRÁFICO DE PARSIMONIA DETALLADO ====================

def plot_parsimony_detail(df, dataset, k_parsimonia, k_best_raw):
    """
    Gráfico detallado que muestra exactamente por qué se eligió ese K:
    - Las 3 métricas individuales
    - El score ponderado
    - La zona de tolerancia del 1%
    - Marcadores claros para K óptimo bruto y K parsimonioso
    """
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    ks_num = list(df.index)
    ks_str = [str(x) for x in ks_num]

    idx_parsimony = ks_num.index(k_parsimony)
    idx_best_raw   = ks_num.index(k_best_raw)

    best_val  = df.loc[k_best_raw, 'Weighted_Score']
    threshold = best_val - PARSIMONY_THRESHOLD

    colors_met = {'BA': '#1f77b4', 'F1': '#ff7f0e', 'PS': '#2ca02c', 'Ponderado': '#9467bd'}

    # ---- Panel 1: Las 3 métricas individuales ----
    ax1 = fig.add_subplot(gs[0, 0])
    for met, col in [('BA', colors_met['BA']),
                     ('F1', colors_met['F1']),
                     ('PS', colors_met['PS'])]:
        ax1.plot(ks_str, df[met], marker='o', color=col, linewidth=2,
                 label=f"{met} (w={METRIC_WEIGHTS[met]:.0%})")

    ax1.axvline(x=idx_parsimony, color='green',  linestyle='-',  alpha=0.6,
                linewidth=2.5, label=f'K elegido ({k_parsimony})')
    ax1.axvline(x=idx_best_raw,   color='red',    linestyle='--', alpha=0.5,
                linewidth=1.5, label=f'K mejor bruto ({k_best_raw})')
    ax1.set_title("Métricas individuales vs K", fontweight='bold')
    ax1.set_xlabel("K (nº características)")
    ax1.set_ylabel("Valor métrica (test)")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_xticklabels(ks_str, rotation=45)

    # ---- Panel 2: Score ponderado + banda de tolerancia ----
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(ks_str, df['Weighted_Score'], marker='D', color=colors_met['Ponderado'],
             linewidth=2.5, markersize=7, label='Score ponderado')

    # Banda de tolerancia (zona verde = dentro del 1%)
    ax2.axhline(y=best_val, color='red',   linestyle='--', alpha=0.6,
                linewidth=1.5, label=f'Máximo ({best_val:.4f})')
    ax2.axhline(y=threshold, color='orange', linestyle='--', alpha=0.6,
                linewidth=1.5, label=f'Umbral parsimonia ({threshold:.4f})')
    ax2.fill_between(range(len(ks_str)), threshold, best_val,
                     alpha=0.15, color='green', label='Zona tolerancia (1%)')

    # Marcador del K elegido
    ax2.scatter([idx_parsimony], [df.loc[k_parsimony, 'Weighted_Score']],
                color='green', s=150, zorder=5, label=f'K elegido ({k_parsimony})')
    ax2.scatter([idx_best_raw],   [df.loc[k_best_raw,   'Weighted_Score']],
                color='red',   s=100, zorder=5, marker='*',
                label=f'K mejor bruto ({k_best_raw})')

    ax2.set_title("Score ponderado + criterio de parsimonia", fontweight='bold')
    ax2.set_xlabel("K (nº características)")
    ax2.set_ylabel("Score ponderado")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(range(len(ks_str)))
    ax2.set_xticklabels(ks_str, rotation=45)

    # ---- Panel 3: Ganancia marginal del score ponderado ----
    ax3 = fig.add_subplot(gs[1, 0])
    weighted_scores = df['Weighted_Score'].values
    deltas     = np.diff(weighted_scores)
    ks_mid     = ks_str[1:]  # un punto menos

    colors_bar = ['#2ca02c' if d > 0 else '#d62728' for d in deltas]
    ax3.bar(ks_mid, deltas, color=colors_bar, alpha=0.75, edgecolor='black', linewidth=0.5)
    ax3.axhline(y=0, color='black', linewidth=0.8)

    # Marcar dónde cae el K parsimonioso
    if idx_parsimony > 0:
        ax3.axvline(x=idx_parsimony - 1, color='green', linestyle='-',
                    alpha=0.6, linewidth=2, label=f'K elegido ({k_parsimony})')
    ax3.set_title("Ganancia marginal del score\n(¿cuánto mejora cada K adicional?)",
                  fontweight='bold')
    ax3.set_xlabel("K (nº características)")
    ax3.set_ylabel("Δ Score ponderado")
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.set_xticklabels(ks_mid, rotation=45)

    # ---- Panel 4: Tabla resumen de las métricas por K ----
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    table_data  = []
    col_labels  = ['K', 'BA', 'F1', 'PS', 'Score\nPond.', 'Elegido']
    cell_colors = []

    for k_idx, k_val in enumerate(ks_num):
        row = df.loc[k_val]
        chosen = 'YES' if k_val == k_parsimony else ''
        table_data.append([
            str(k_val),
            f"{row['BA']:.4f}",
            f"{row['F1']:.4f}",
            f"{row['PS']:.4f}",
            f"{row['Weighted_Score']:.4f}",
            chosen
        ])
        if k_val == k_parsimony:
            cell_colors.append(['#d4f1be'] * 6)
        elif k_val == k_best_raw:
            cell_colors.append(['#ffd6d6'] * 6)
        else:
            cell_colors.append(['white'] * 6)

    table = ax4.table(cellText=table_data, colLabels=col_labels,
                      cellLoc='center', loc='center',
                      cellColours=cell_colors)
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.1, 1.4)
    ax4.set_title("Tabla completa de métricas por K\n(verde=elegido, rojo=mejor bruto)",
                  fontweight='bold', pad=10)

    fig.suptitle(f"Análisis de Parsimonia — {dataset.upper()}\n"
                  f"Pesos: BA={METRIC_WEIGHTS['BA']:.0%}, "
                  f"F1={METRIC_WEIGHTS['F1']:.0%}, PS={METRIC_WEIGHTS['PS']:.0%}",
                  fontsize=13, fontweight='bold')

    out_path = os.path.join(ANALYSIS_SUBDIR, "plots",
                            f"parsimony_detail_{dataset}.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     Gráfico detallado guardado: {out_path}")

# ==================== VISUALIZACIÓN: TABLA RESUMEN GLOBAL ====================

def save_summary_table(summary_rows):
    """Guarda la tabla resumen en CSV, en texto detallado y el informe ejecutivo."""

    df_summary = pd.DataFrame(summary_rows)

    # 1. CSV
    csv_path = os.path.join(ANALYSIS_SUBDIR, "parsimony_summary.csv")
    df_summary.to_csv(csv_path, index=False)

    # 2. Texto detallado de parsimonia
    txt_path = os.path.join(ANALYSIS_SUBDIR, "parsimony_summary.txt")
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("INFORME DE SELECCIÓN POR PARSIMONIA\n")
        f.write(f"Pesos: BA={METRIC_WEIGHTS['BA']:.0%}  "
                f"F1={METRIC_WEIGHTS['F1']:.0%}  "
                f"PS={METRIC_WEIGHTS['PS']:.0%}\n")
        f.write(f"Umbral tolerancia: {PARSIMONY_THRESHOLD*100:.1f}%\n")
        f.write("=" * 80 + "\n\n")

        for row in summary_rows:
            f.write(f"DATASET: {row['Dataset'].upper()}\n")
            f.write(f"  K elegido (parsimonia): {row['Chosen_K']}\n")
            f.write(f"  K mejor bruto:          {row['Raw_Best_K']}\n")
            f.write(f"  Genes seleccionados:    {row['Num_Genes']}\n")
            f.write(f"  Metodología ganadora:   {row['Methodology']}\n")
            f.write(f"  BA  (test):             {row['BA']:.4f}\n")
            f.write(f"  F1  (test):             {row['F1']:.4f}\n")
            f.write(f"  PS  (test):             {row['PS']:.4f}\n")
            f.write(f"  Score ponderado:        {row['Weighted_Score']:.4f}\n")
            f.write(f"  Score mejor bruto:      {row['Raw_Best_Score']:.4f}\n")
            f.write(f"  Diferencia:             "
                    f"{row['Raw_Best_Score'] - row['Weighted_Score']:.4f} "
                    f"({'dentro del umbral' if abs(row['Raw_Best_Score'] - row['Weighted_Score']) <= PARSIMONY_THRESHOLD else 'SUPERA UMBRAL'})\n")
            f.write("-" * 60 + "\n\n")

    # 3. INFORME EJECUTIVO (Pedido por el usuario)
    exec_path = os.path.join(final_output_root, "best_winner_summary.txt")
    with open(exec_path, 'w', encoding='utf-8') as f:
        f.write("INFORME EJECUTIVO DE RESULTADOS FINALES\n")
        f.write("=" * 60 + "\n\n")
        for row in summary_rows:
            f.write(f"DATASET: {row['Dataset'].upper()}\n")
            f.write(f"  K Final: {row['Chosen_K']}\n")
            f.write(f"  Metodología: {row['Methodology']}\n")
            f.write(f"  Scores: BA={row['BA']:.4f}, F1={row['F1']:.4f}, PS={row['PS']:.4f}\n")
            f.write("-" * 40 + "\n")

    # Tabla visual en consola
    print("\n" + "=" * 80)
    print(" TABLA RESUMEN FINAL")
    print("=" * 80)
    print(df_summary.to_string(index=False))

    print(f"\n    -> CSV guardado:        {csv_path}")
    print(f"    -> Texto detallado:     {txt_path}")
    print(f"    -> INFORME EJECUTIVO:   {exec_path}")

# ==================== PROCESO PRINCIPAL ====================

setup_directories()

print("=" * 80)
print(" CONSOLIDACIÓN FINAL: PARSIMONIA, DATOS Y GRÁFICOS")
print(f"   Pesos: BA={METRIC_WEIGHTS['BA']:.0%}, "
      f"F1={METRIC_WEIGHTS['F1']:.0%}, PS={METRIC_WEIGHTS['PS']:.0%}")
print(f"   Umbral parsimonia: {PARSIMONY_THRESHOLD*100:.1f}%")
print("=" * 80)

all_winners_data = []
summary_rows     = []

for dataset in datasets:
    print(f"\n Analizando {dataset.upper()}...")
    rows = []

    for k in k_values:
        py_path = os.path.join(results_root, dataset, f"k_{k}",
                               "best_analysis",
                               f"{dataset}_best_dataset_k{k}.py")
        info = parse_best_dataset_py(py_path)

        if info:
            score = weighted_score(info['BA_Test'], info['F1_Test'], info['PS_Test'])
            rows.append({
                'K':         k,
                'Methodology': info['Ganador'],
                'BA':        info['BA_Test'],
                'F1':        info['F1_Test'],
                'PS':        info['PS_Test'],
                'Weighted_Score': score
            })
            all_winners_data.append({
                'Dataset': dataset, 'K': k, 'Metodo': info['Ganador']
            })

    if not rows:
        print(f"    No se encontraron datos para {dataset}")
        continue

    df = pd.DataFrame(rows).set_index('K')

    # --- Selección por parsimonia ---
    k_best_raw    = df['Weighted_Score'].idxmax()
    best_val      = df.loc[k_best_raw, 'Weighted_Score']
    threshold     = best_val - PARSIMONY_THRESHOLD
    k_parsimony  = df[df['Weighted_Score'] >= threshold].index.min()
    row_winner    = df.loc[k_parsimony]

    print(f"    K mejor bruto: {k_best_raw} (score={best_val:.4f})")
    print(f"    K elegido por parsimonia: {k_parsimony} "
          f"(score={row_winner['Weighted_Score']:.4f}, "
          f"diff={best_val - row_winner['Weighted_Score']:.4f})")

    # --- Guardar tabla de métricas por K como CSV ---
    df_export = df.reset_index().rename(columns={
        'K':         'K',
        'BA':        'BA_Test',
        'F1':        'F1_Test',
        'PS':        'PS_Test'
    })
    df_export['Chosen'] = df_export['K'].apply(
        lambda k: 'YES' if k == k_parsimony else '')
    df_export['Raw_Best_K'] = df_export['K'].apply(
        lambda k: 'YES' if k == k_best_raw else '')

    csv_k_path = os.path.join(ANALYSIS_SUBDIR, f"metrics_per_k_{dataset}.csv")
    df_export.to_csv(csv_k_path, index=False)
    print(f"    Tabla métricas por K guardada: {csv_k_path}")

    # --- Gráfico detallado de parsimonia ---
    plot_parsimony_detail(df, dataset, k_parsimony, k_best_raw)

    # --- Rescate de datos y genes ---
    winner_folder = os.path.join(final_output_root,
                                 f"{dataset}_k{k_parsimony}_winner")
    os.makedirs(winner_folder, exist_ok=True)

    src_k_folder = os.path.join(results_root, dataset, f"k_{k_parsimony}")
    train_name, test_name = get_dataset_filename_pattern(
        dataset, k_parsimony, row_winner['Methodology'])

    for orig, dest in [(train_name, "train_dataset.csv"),
                       (test_name,  "test_dataset.csv")]:
        src_path = os.path.join(src_k_folder, orig)
        if os.path.exists(src_path):
            shutil.copy(src_path, os.path.join(winner_folder, dest))
        else:
            print(f"      No encontrado: {src_path}")

    # Extracción de genes
    n_genes = 0
    train_dest = os.path.join(winner_folder, "train_dataset.csv")
    if os.path.exists(train_dest):
        try:
            df_cols = pd.read_csv(train_dest, nrows=0)
            genes   = [c for c in df_cols.columns if c.lower() != 'target']
            n_genes = len(genes)
            with open(os.path.join(winner_folder, "selected_genes.txt"), "w") as fg:
                fg.write("\n".join(genes))
            print(f"    {n_genes} genes exportados.")
        except Exception as e:
            print(f"    Error al extraer genes: {e}")

    # Acumular para tabla resumen
    summary_rows.append({
        'Dataset':           dataset,
        'Chosen_K':         k_parsimony,
        'Raw_Best_K':     k_best_raw,
        'Num_Genes':           n_genes,
        'Methodology':       row_winner['Methodology'],
        'BA':                round(row_winner['BA'],         4),
        'F1':                round(row_winner['F1'],         4),
        'PS':                round(row_winner['PS'],         4),
        'Weighted_Score':   round(row_winner['Weighted_Score'],  4),
        'Raw_Best_Score': round(best_val,                 4)
    })

# ==================== TABLA Y TEXTO RESUMEN GLOBAL ====================

if summary_rows:
    save_summary_table(summary_rows)

# ==================== GRÁFICO GLOBAL: DISTRIBUCIÓN DE MÉTODOS ====================

if all_winners_data:
    plt.figure(figsize=(12, 8))
    sns.countplot(data=pd.DataFrame(all_winners_data),
                  y='Metodo', hue='Dataset', palette='viridis')
    plt.title("Frecuencia de Metodologías Ganadoras por Dataset",
              fontsize=13, fontweight='bold')
    plt.xlabel("Veces que fue mejor en un valor de K")
    plt.ylabel("Algoritmo / Metodología")
    plt.tight_layout()
    out_path = os.path.join(ANALYSIS_SUBDIR, "plots", "methods_distribution.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"\n    Distribución de métodos guardada: {out_path}")

print("\n" + "=" * 80)
print(f" CONSOLIDACIÓN FINALIZADA.")
print(f"    Resultados principales: '{final_output_root}/'")
print(f"    Análisis de parsimonia: '{ANALYSIS_SUBDIR}/'")
print("=" * 80)