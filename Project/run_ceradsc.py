import subprocess
import os
import time
import shutil
import glob

# --- CONFIGURACIÓN PARA CERADSC ---
DATASET_NAME = "ceradsc"
k_values = [25, 50, 75, 100, 500, 1000, 1500, 2500]

# Obtenemos la ruta absoluta de la carpeta donde está este script (Project/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SCRIPTS_DIR = os.path.join(SCRIPT_DIR, "scripts")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results_ml_classification")

script_baseline = "1.1-imbalanced-classification-readcsv.py"
baseline_dummy_k = 25 
iterative_scripts = [
    "1.2-imbalanced-classification-readcsv-FR.py",
    "1.3-Re-Feature-imbalanced-classification-readcsv-FR.py"
]
judge_script = "4-Statistical_tests.py"

# --- FUNCIONES ---
def run_command(script_name, dataset, k, input_dir, output_dir):
    script_path = os.path.join(SCRIPTS_DIR, script_name)
    if not os.path.exists(script_path):
        print(f" Error: No encuentro {script_path}")
        return False
    cmd = ["python", script_path, "--dataset", dataset, "--k", str(k), "--input_dir", input_dir, "--output_root", output_dir]
    try:
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f" FALLO en {script_name} ({dataset} k={k}): {e}")
        return False

def distribute_baseline(dataset, source_k, target_k):
    src_dir = os.path.join(RESULTS_DIR, dataset, f"k_{source_k}")
    dst_dir = os.path.join(RESULTS_DIR, dataset, f"k_{target_k}")
    if not os.path.exists(src_dir): return
    os.makedirs(dst_dir, exist_ok=True)
    for file_path in glob.glob(os.path.join(src_dir, "*.csv")):
        filename = os.path.basename(file_path)
        if f"k{source_k}" in filename:
            new_name = filename.replace(f"k{source_k}", f"k{target_k}")
            shutil.copy(file_path, os.path.join(dst_dir, new_name))


# --- MAIN ---
if __name__ == "__main__":
    print(f" INICIANDO AUTOMATIZACIÓN PARA: {DATASET_NAME.upper()}")
    
    # 1. Baseline (Una sola vez)
    print(" Ejecutando Baseline 1.1...")
    if run_command(script_baseline, DATASET_NAME, baseline_dummy_k, DATA_DIR, RESULTS_DIR):
        # 2. Bucle K
        for k in k_values:
            print(f"\n Procesando K={k}...")
            os.makedirs(os.path.join(RESULTS_DIR, DATASET_NAME, f"k_{k}"), exist_ok=True)
            
            if k != baseline_dummy_k:
                distribute_baseline(DATASET_NAME, baseline_dummy_k, k)
            
            for script in iterative_scripts:
                print(f"    {script}")
                run_command(script, DATASET_NAME, k, DATA_DIR, RESULTS_DIR)
            
            print(f"    Juez...")
            run_command(judge_script, DATASET_NAME, k, RESULTS_DIR, RESULTS_DIR)
            
        print(f" FINALIZADO {DATASET_NAME}")