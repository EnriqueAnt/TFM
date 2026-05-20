import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import seaborn as sns
import itertools

# ==================== CONFIGURACIÓN ====================
# Obtenemos la ruta absoluta de la carpeta donde está este script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

RESULTS_DIR = os.path.join(SCRIPT_DIR, "best_k_selection")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "biological_analysis")
datasets_target = ["braaksc", "ceradsc", "cogdx"]

# Crear estructura de salida
if os.path.exists(OUTPUT_DIR):
    import shutil
    shutil.rmtree(OUTPUT_DIR)
os.makedirs(OUTPUT_DIR)
os.makedirs(os.path.join(OUTPUT_DIR, "plots"))
os.makedirs(os.path.join(OUTPUT_DIR, "gene_lists"))

# ==================== FUNCIONES AUXILIARES ====================

def find_gene_file(dataset_name):
    """Busca el archivo selected_genes.txt en la carpeta del ganador"""
    for folder in os.listdir(RESULTS_DIR):
        if dataset_name in folder and os.path.isdir(os.path.join(RESULTS_DIR, folder)):
            file_path = os.path.join(RESULTS_DIR, folder, "selected_genes.txt")
            if os.path.exists(file_path):
                return file_path
    return None

def load_genes(dataset_name):
    path = find_gene_file(dataset_name)
    if path:
        with open(path, 'r') as f:
            genes = set(line.strip() for line in f if line.strip())
        print(f"  {dataset_name.upper()}: Cargados {len(genes)} genes")
        return genes
    else:
        print(f"  {dataset_name.upper()}: No se encontró lista de genes")
        return set()

def save_list(name, gene_set, description=""):
    """Guarda una lista de genes en un archivo y la imprime en consola"""
    filename = os.path.join(OUTPUT_DIR, "gene_lists", f"{name}.txt")
    genes_sorted = sorted(list(gene_set))
    
    with open(filename, "w") as f:
        f.write("\n".join(genes_sorted))
    
    print(f"\n  {description.upper()} ({len(genes_sorted)} genes):")
    if genes_sorted:
        # Imprimimos primeros 50 para no saturar, pero el archivo los tiene todos
        preview = ", ".join(genes_sorted[:50])
        dots = "..." if len(genes_sorted) > 50 else ""
        print(f"   -> {preview}{dots}")
    else:
        print("   -> Ninguno")
    print(f"   (Guardado en: {filename})")

def jaccard_index(set_a, set_b):
    if not set_a and not set_b: return 0.0
    return len(set_a.intersection(set_b)) / len(set_a.union(set_b))

def overlap_ratio(set_a, set_b):
    if not set_a:
        return 0.0
    return len(set_a.intersection(set_b)) / len(set_a)
    
# ==================== PROCESO PRINCIPAL ====================

print("=" * 60)
print("  ANÁLISIS DE INTERSECCIONES BIOLÓGICAS (COMPLETO)")
print("=" * 60)

# 1. Cargar Genes
gene_sets = {}
for ds in datasets_target:
    gene_sets[ds] = load_genes(ds)

if len(gene_sets) < 3:
    print("  Faltan datasets. Se requieren los 3 para este análisis.")
    exit()

# ---------------------------------------------------------
# A. ANÁLISIS DE SIMILITUD (JACCARD)
# ---------------------------------------------------------
print("\n" + "-"*60)
print("A. MATRIZ DE SIMILITUD")
print("-" * 60)
jaccard_data = []
pairs = list(itertools.combinations(datasets_target, 2))

print(f"{'Comparación':<20} | {'Jaccard':<10} | {'A->B':<10} | {'B->A':<10} | {'Comunes':<10}")
for d1, d2 in pairs:
    jacc = jaccard_index(gene_sets[d1], gene_sets[d2])
    comunes = len(gene_sets[d1] & gene_sets[d2])
    overlap_ab = overlap_ratio(gene_sets[d1], gene_sets[d2])
    overlap_ba = overlap_ratio(gene_sets[d2], gene_sets[d1])
    print(f"{d1.upper()} vs {d2.upper():<10} | {jacc:.4f}     | {overlap_ab:.4f}   | {overlap_ba:.4f}   | {comunes}")
    

    jaccard_data.append({
        'Dataset_A': d1.upper(),
        'Dataset_B': d2.upper(),
        'Jaccard_Index': jacc,
        'Ratio_Inclusion_A_en_B': overlap_ab,
        'Ratio_Inclusion_B_en_A': overlap_ba,
        'Common_Genes_Count': comunes
    })

# Guardar el archivo 
pd.DataFrame(jaccard_data).to_csv(os.path.join(OUTPUT_DIR, "intersection_metrics_detailed.csv"), index=False)


# Matriz completa de Jaccard
matrix = pd.DataFrame(index=datasets_target, columns=datasets_target)

for d1 in datasets_target:
    for d2 in datasets_target:
        matrix.loc[d1, d2] = jaccard_index(gene_sets[d1], gene_sets[d2])

matrix.to_csv(os.path.join(OUTPUT_DIR, "jaccard_matrix.csv"))

plt.figure(figsize=(6,5))
sns.heatmap(matrix.astype(float), annot=True, fmt=".2f")
plt.title("Jaccard Similarity Matrix")
plt.savefig(os.path.join(OUTPUT_DIR, "plots", "heatmap_jaccard.png"), dpi=300)
plt.close()

# ---------------------------------------------------------
# B. GENES MAESTROS (CORE DRIVERS)
# ---------------------------------------------------------
print("\n" + "-"*60)
print("B. GENES MAESTROS (INTERSECCIÓN TRIPLE)")
print("-" * 60)
core_genes = gene_sets['braaksc'] & gene_sets['ceradsc'] & gene_sets['cogdx']
save_list("01_genes_core_masters", core_genes, "Genes Maestros (Braak + Cerad + Cogdx)")

# ---------------------------------------------------------
# C. GENES EXCLUSIVOS (BIOMARCADORES ESPECÍFICOS)
# ---------------------------------------------------------
print("\n" + "-"*60)
print("C. GENES EXCLUSIVOS (UNIQUE)")
print("-" * 60)

braak_only = gene_sets['braaksc'] - gene_sets['ceradsc'] - gene_sets['cogdx']
save_list("02_genes_unique_braak", braak_only, "Exclusivos Braak (Solo Tau)")

cerad_only = gene_sets['ceradsc'] - gene_sets['braaksc'] - gene_sets['cogdx']
save_list("03_genes_unique_cerad", cerad_only, "Exclusivos Cerad (Solo Amiloide)")

cogdx_only = gene_sets['cogdx'] - gene_sets['braaksc'] - gene_sets['ceradsc']
save_list("04_genes_unique_cogdx", cogdx_only, "Exclusivos Cogdx (Clínica sin Patología Clásica)")

# ---------------------------------------------------------
# D. CRUCES CLÍNICOS ESPECÍFICOS (PARES VS UNO)
# ---------------------------------------------------------
print("\n" + "-"*60)
print("D. DISOCIACIONES CLÍNICAS (2 vs 1)")
print("-" * 60)

# 1. RESILIENCIA: Tienen Braak y Cerad, pero NO Cogdx
# (Cerebro dañado físicamente, pero paciente sano cognitivamente)
pathology_no_symptoms = (gene_sets['braaksc'] & gene_sets['ceradsc']) - gene_sets['cogdx']
save_list("05_genes_pathology_resilience", pathology_no_symptoms, "Resiliencia (Braak + Cerad - Cogdx)")

# 2. TAU-DEPENDIENTE: Tienen Braak y Cogdx, pero NO Cerad
# (Síntomas explicados por Tau, sin influencia de Amiloide)
tau_cognitive = (gene_sets['braaksc'] & gene_sets['cogdx']) - gene_sets['ceradsc']
save_list("06_genes_tau_cognitive", tau_cognitive, "Tau-Cognitivo (Braak + Cogdx - Cerad)")

# 3. AMILOIDE-DEPENDIENTE: Tienen Cerad y Cogdx, pero NO Braak
# (Síntomas explicados por Amiloide, sin influencia de Tau)
amyloid_cognitive = (gene_sets['ceradsc'] & gene_sets['cogdx']) - gene_sets['braaksc']
save_list("07_genes_amyloid_cognitive", amyloid_cognitive, "Amiloide-Cognitivo (Cerad + Cogdx - Braak)")

# ---------------------------------------------------------
# E. GRÁFICOS
# ---------------------------------------------------------
print("\n  Generando Diagramas de Venn...")

# Venn Triple
plt.figure(figsize=(10, 10))
try:
    venn3([gene_sets['braaksc'], gene_sets['ceradsc'], gene_sets['cogdx']], 
          set_labels=('Braak', 'Cerad', 'Cogdx'))
    plt.title("Mapa Global de Intersecciones")
    plt.savefig(os.path.join(OUTPUT_DIR, "plots", "venn_triple.png"), dpi=300)
except Exception as e:
    print(f"No se pudo generar Venn Triple: {e}")
plt.close()

# Venns Pares
for d1, d2 in pairs:
    plt.figure(figsize=(8, 8))
    try:
        venn2([gene_sets[d1], gene_sets[d2]], set_labels=(d1.upper(), d2.upper()))
        plt.title(f"Intersección: {d1.upper()} vs {d2.upper()}")
        plt.savefig(os.path.join(OUTPUT_DIR, "plots", f"venn_{d1}_{d2}.png"), dpi=300)
    except: pass
    plt.close()

print(f"\n  ANÁLISIS COMPLETADO. Archivos guardados en '{OUTPUT_DIR}/gene_lists'")