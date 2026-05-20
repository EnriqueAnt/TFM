import os
import pandas as pd
import glob
import re
import json
import io

def get_selected_features(base_dir, dataset, k_val, methodology):
    """
    Recupera la lista de genes seleccionados para una metodología y dataset dados.
    Prioridad: metadatos JSON -> archivos de entrenamiento CSV.
    """
    metadata_dir = os.path.join(base_dir, dataset, k_val, "runs_metadata")
    
    # Limpiar nombre de metodología para mapeo JSON
    clean_meth = methodology.replace("_Resampl", "").replace("_FR", "").replace("_Baseline", "")
    json_filename = f"genes-{dataset}_{clean_meth}_run1.json"
    json_path = os.path.join(metadata_dir, json_filename)
    
    if os.path.exists(json_path):
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                if "selected_features" in data:
                    return set(data["selected_features"])
        except Exception:
            pass

    # Fallback to CSV if JSON is missing (e.g., R1 or legacy runs)
    k_num = k_val.split('_')[1]
    branch = "FR"
    if "_Resampl" in methodology or "R3_Resampl" in methodology:
        branch = "resampling_FR"
    elif "_Baseline" in methodology or "R1_Baseline" in methodology:
        branch = "Baseline"
    
    if branch == "Baseline":
        filename = f"genes-{dataset}_train.csv"
    else:
        filename = f"genes-{dataset}-{branch}-{k_num}_train.csv"
        
    file_path = os.path.join(base_dir, dataset, k_val, filename)
    
    if os.path.exists(file_path):
        try:
            cols = pd.read_csv(file_path, nrows=0).columns.tolist()
            if "target" in cols: cols.remove("target")
            return set(cols)
        except Exception:
            return None
    return None

def run_methodology_audit():
    """
    Función de auditoría principal para comparar ganadores originales vs ganadores ponderados personalizados.
    Genera reportes y resúmenes CSV.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    base_dir = os.path.join(project_root, "results_ml_classification")
    output_dir = os.path.join(script_dir, "robustness_audit_results")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # --- Cargar ganadores globales del reporte final ---
    global_winners = {}
    report_path = os.path.join(project_root, "best_k_selection", "best_winner_summary.txt")
    if os.path.exists(report_path):
        with open(report_path, "r", encoding="utf-8") as f:
            content = f.read()
            blocks = re.split(r"-{10,}", content)
            for block in blocks:
                ds_match = re.search(r"DATASET:\s*(\w+)", block, re.I)
                k_match = re.search(r"K Final:\s*(\d+)", block)
                meth_match = re.search(r"Metodología:\s*([^\n\r]+)", block)
                if ds_match and k_match and meth_match:
                    global_winners[ds_match.group(1).lower()] = {
                        "k": f"k_{k_match.group(1)}",
                        "meth": meth_match.group(1).strip()
                    }

    if not os.path.exists(base_dir):
        print(f"Error: Results directory not found at {base_dir}")
        return

    # Find Signed Rank summary files
    summary_files = glob.glob(os.path.join(base_dir, "*", "k_*", "best_analysis", "summary_signed_rank.csv"))
    if not summary_files:
        print("No results found in the project tree.")
        return

    # Ordenar archivos por dataset y valor k
    def sort_key(path):
        p = path.split(os.sep)
        return (p[-4], int(p[-3].split('_')[1]))
    summary_files.sort(key=sort_key)

    sr_results = []
    mean_results = []
    detailed_comparisons = []
    methodology_mismatches = {} # (dataset, k_val, table_type) -> (orig, pond)
    local_winners_map = {} 

    for sr_path in summary_files:
        parts = sr_path.split(os.sep)
        dataset = parts[-4]
        k_val = parts[-3]
        mean_path = sr_path.replace("summary_signed_rank.csv", "summary_mean.csv")

        if not os.path.exists(mean_path): continue

        df_sr = pd.read_csv(sr_path, index_col=0)
        df_mean = pd.read_csv(mean_path, index_col=0)

        # Custom Score Calculation (50% BA, 35% F1, 15% PS)
        w_ba, w_f1, w_ps = 0.5, 0.35, 0.15
        df_sr['Custom_Score'] = (df_sr['BA'] * w_ba) + (df_sr['F1'] * w_f1) + (df_sr['PS'] * w_ps)
        df_mean['Custom_Score'] = (df_mean['BA'] * w_ba) + (df_mean['F1'] * w_f1) + (df_mean['PS'] * w_ps)

        # Identify best methodologies
        best_sr_total = df_sr['Total'].idxmax()
        best_sr_custom = df_sr['Custom_Score'].idxmax()
        best_mean_total = df_mean['Total'].idxmax()
        best_mean_custom = df_mean['Custom_Score'].idxmax()

        # Store Table 1 & 2 data
        sr_results.append({
            "Dataset": dataset, "K": k_val, 
            "Best_Original": best_sr_total, "Best_Weighted": best_sr_custom
        })
        mean_results.append({
            "Dataset": dataset, "K": k_val, 
            "Best_Original": best_mean_total, "Best_Weighted": best_mean_custom
        })
        
        # Capturar discrepancias para la Tabla 3 detallada y análisis de genes
        if best_sr_total != best_sr_custom:
            methodology_mismatches[(dataset, k_val, "SR")] = (best_sr_total, best_sr_custom)
            for meth, label in [(best_sr_total, "Orig"), (best_sr_custom, "Pond")]:
                detailed_comparisons.append({
                    "Table": "SR", "Dataset": dataset, "K": k_val, "Type": label, "Methodology": meth,
                    "BA": df_sr.loc[meth, "BA"], "F1": df_sr.loc[meth, "F1"], "PS": df_sr.loc[meth, "PS"],
                    "Total": df_sr.loc[meth, "Total"], "Custom": df_sr.loc[meth, "Custom_Score"]
                })

        if best_mean_total != best_mean_custom:
            methodology_mismatches[(dataset, k_val, "Mean")] = (best_mean_total, best_mean_custom)
            for meth, label in [(best_mean_total, "Orig"), (best_mean_custom, "Pond")]:
                detailed_comparisons.append({
                    "Table": "Mean", "Dataset": dataset, "K": k_val, "Type": label, "Methodology": meth,
                    "BA": df_mean.loc[meth, "BA"], "F1": df_mean.loc[meth, "F1"], "PS": df_mean.loc[meth, "PS"],
                    "Total": df_mean.loc[meth, "Total"], "Custom": df_mean.loc[meth, "Custom_Score"]
                })

        local_winners_map[(dataset, k_val)] = {
            "SR_T": best_sr_total, "SR_P": best_sr_custom,
            "M_T": best_mean_total, "M_P": best_mean_custom
        }

    # --- GUARDAR ARCHIVOS CSV ---
    pd.DataFrame(sr_results).to_csv(os.path.join(output_dir, "table1_sr_performance.csv"), index=False)
    pd.DataFrame(mean_results).to_csv(os.path.join(output_dir, "table2_mean_performance.csv"), index=False)
    if detailed_comparisons:
        pd.DataFrame(detailed_comparisons).to_csv(os.path.join(output_dir, "table3_detailed_divergence.csv"), index=False)

    # --- GENERAR REPORTE COMPLETO EN TEXTO ---
    report = io.StringIO()
    
    def rprint(text=""):
        print(text)
        report.write(text + "\n")

    rprint("=" * 130)
    rprint(f"{'REPORT: METHODOLOGY ROBUSTNESS AND DIVERGENCE ANALYSIS':^130}")
    rprint("=" * 130)

    # TABLA 1
    rprint("\n" + "=" * 130)
    rprint(f"{'TABLE 1: SIGNED RANK (SR) ANALYSIS':^130}")
    rprint("=" * 130)
    rprint(f"{'Dataset':<15} | {'K':<7} | {'Best Methodology (Total Score)':<45} | {'Best Methodology (Weighted Score)'}")
    rprint("-" * 130)
    for r in sr_results:
        rprint(f"{r['Dataset']:<15} | {r['K']:<7} | {r['Best_Original']:<45} | {r['Best_Weighted']}")
    
    m_sr = [f"{d} ({k})" for (d, k, t), p in methodology_mismatches.items() if t == "SR"]
    rprint("-" * 130)
    rprint(f"Mismatches detected in SR: {', '.join(m_sr) if m_sr else 'None'}")

    # TABLA 2
    rprint("\n\n" + "=" * 130)
    rprint(f"{'TABLE 2: MEAN PERFORMANCE ANALYSIS':^130}")
    rprint("=" * 130)
    rprint(f"{'Dataset':<15} | {'K':<7} | {'Best Methodology (Total Score)':<45} | {'Best Methodology (Weighted Score)'}")
    rprint("-" * 130)
    for r in mean_results:
        rprint(f"{r['Dataset']:<15} | {r['K']:<7} | {r['Best_Original']:<45} | {r['Best_Weighted']}")
    
    m_mean = [f"{d} ({k})" for (d, k, t), p in methodology_mismatches.items() if t == "Mean"]
    rprint("-" * 130)
    rprint(f"Mismatches detected in Mean: {', '.join(m_mean) if m_mean else 'None'}")

    # TABLA 3
    if detailed_comparisons:
        rprint("\n\n" + "=" * 165)
        rprint(f"{'TABLE 3: DETAILED PERFORMANCE METRICS FOR MISMATCHES':^165}")
        rprint("=" * 165)
        rprint(f"{'Table':<5} | {'Dataset':<10} | {'K':<10} | {'Type':<5} | {'BA':<4} | {'F1':<4} | {'PS':<4} | "
               f"{'Total':<6} | {'Wght.':<6} | {'W_SR_T':<6} | {'W_SR_P':<6} | {'W_M_T':<6} | {'W_M_P':<6} | {'Methodology'}")
        rprint("-" * 165)
        
        last_group = None
        for d in detailed_comparisons:
            ds, k_val = d['Dataset'], d['K']
            if last_group and (ds, k_val) != last_group: rprint("-" * 165)
            
            winners = local_winners_map.get((ds, k_val), {})
            meth = d['Methodology']
            
            def get_mark(meth, win_meth, col_type):
                if d['Type'] == col_type: return "-"
                return "YES" if meth == win_meth else "no"

            w_sr_t = get_mark(meth, winners.get("SR_T"), "Orig" if d['Table'] == 'SR' else None)
            w_sr_p = get_mark(meth, winners.get("SR_P"), "Pond" if d['Table'] == 'SR' else None)
            w_m_t = get_mark(meth, winners.get("M_T"), "Orig" if d['Table'] == 'Mean' else None)
            w_m_p = get_mark(meth, winners.get("M_P"), "Pond" if d['Table'] == 'Mean' else None)

            rprint(f"{d['Table']:<5} | {ds:<10} | {k_val:<10} | {d['Type']:<5} | "
                   f"{d['BA']:<4.0f} | {d['F1']:<4.0f} | {d['PS']:<4.0f} | "
                   f"{d['Total']:<6.0f} | {d['Custom']:<6.1f} | "
                   f"{w_sr_t:<6} | {w_sr_p:<6} | {w_m_t:<6} | {w_m_p:<6} | {meth}")
            last_group = (ds, k_val)

    # ANÁLISIS DE GENES
    if methodology_mismatches:
        rprint("\n\n" + "=" * 130)
        rprint(f"{'GENE SELECTION STABILITY ANALYSIS (ORIGINAL VS WEIGHTED)':^130}")
        rprint("=" * 130)
        
        global_gene_sets = {}
        for ds_name, info in global_winners.items():
            genes = get_selected_features(base_dir, ds_name.upper(), info['k'], info['meth'])
            if genes: global_gene_sets[ds_name] = genes
        
        all_sets = list(global_gene_sets.values())
        global_intersection = set.intersection(*all_sets) if len(all_sets) == 3 else set()
        global_union = set.union(*all_sets) if all_sets else set()

        for (ds, k, table_type), (m1, m2) in methodology_mismatches.items():
            genes1 = get_selected_features(base_dir, ds, k, m1)
            genes2 = get_selected_features(base_dir, ds, k, m2)
            
            rprint(f"[{table_type}] Dataset: {ds:<10} | K: {k:<6}")
            rprint(f"  - Original Method: {m1}")
            rprint(f"  - Weighted Method: {m2}")
            
            if genes1 is None or genes2 is None:
                rprint("  - WARNING: Could not find feature lists for one or both methods.")
            elif genes1 == genes2:
                rprint("  - RESULT: Gene selection is IDENTICAL across both methodologies.")
            else:
                g_in = genes2 - genes1
                rprint(f"  - RESULT: Gene selection CHANGED. {len(g_in)} genes substituted (Total K: {len(genes1)}).")
            
            # Análisis de intersección y unión
            rprint(f"  - MULTI-DATASET INTEGRATION ANALYSIS:")
            
            # vs Ganador Propio
            own_info = global_winners.get(ds.lower())
            if own_info:
                genes_own = global_gene_sets.get(ds.lower())
                if genes_own:
                    rprint(f"    * vs Global Winner ({own_info['k']} | {own_info['meth']}):")
                    for name, genes_curr in [("Original", genes1), ("Weighted", genes2)]:
                        if genes_curr == genes_own: res = "IDENTICAL"
                        elif len(genes_curr) < len(genes_own):
                            res = "Subset (Contained)" if genes_curr.issubset(genes_own) else "Different"
                        else:
                            res = "Superset (Includes Winner)" if genes_own.issubset(genes_curr) else "Different"
                        rprint(f"      - {name:<9} Methodology: {res}")

            # vs Intersección
            if global_intersection:
                rprint(f"    * vs GLOBAL INTERSECTION (Common Nucleus | {len(global_intersection)} genes):")
                for name, genes_curr in [("Original", genes1), ("Weighted", genes2)]:
                    res = "CONTAINS CONSENSUS (Full)" if global_intersection.issubset(genes_curr) else f"PARTIAL (Missing {len(global_intersection - genes_curr)} core genes)"
                    rprint(f"      - {name:<9} Methodology: {res}")

            # vs Unión
            if global_union:
                rprint(f"    * vs GLOBAL UNION (Project Frontier | {len(global_union)} genes):")
                for name, genes_curr in [("Original", genes1), ("Weighted", genes2)]:
                    res = "CONSISTENT (All genes within Union)" if genes_curr.issubset(global_union) else f"EXTENDED ({len(genes_curr - global_union)} new genes outside Union)"
                    rprint(f"      - {name:<9} Methodology: {res}")
            rprint("-" * 130)

    # ANÁLISIS DE DIVERGENCIA FINAL (DEDUPLICADO)
    rprint("\n" + "=" * 130)
    rprint(f"{'DIVERGENCE ANALYSIS: WEIGHTED METHODS VS GLOBAL STANDARDS':^130}")
    rprint("=" * 130)
    
    unique_ponds = set((ds, k, m) for (ds, k, t), (o, m) in methodology_mismatches.items())
    for ds_name, k_val, m_pond in sorted(list(unique_ponds)):
        genes_pond = get_selected_features(base_dir, ds_name, k_val, m_pond)
        if not genes_pond: continue
        
        win_info = global_winners.get(ds_name.lower())
        genes_win = get_selected_features(base_dir, ds_name.upper(), win_info['k'], win_info['meth']) if win_info else None
        
        rprint(f"Dataset: {ds_name:<10} | K: {k_val:<6} | Weighted Method: {m_pond}")
        
        if genes_win:
            plus, minus = sorted(list(genes_pond - genes_win)), sorted(list(genes_win - genes_pond))
            rprint(f"  - vs GLOBAL WINNER ({win_info['k']} | {win_info['meth']}):")
            if plus:
                rprint(f"    * ADDITIONAL Genes (New features) [{len(plus)}]:")
                for i in range(0, len(plus), 10): rprint(f"      {', '.join(plus[i:i+10])}")
            if minus:
                rprint(f"    * MISSING Genes (Lost features) [{len(minus)}]:")
                for i in range(0, len(minus), 10): rprint(f"      {', '.join(minus[i:i+10])}")
        
        if global_union:
            u_plus = sorted(list(genes_pond - global_union))
            u_minus = sorted(list(global_union - genes_pond))
            rprint(f"  - vs GLOBAL UNION (All Winning Features):")
            if u_plus:
                rprint(f"    * EXTERNAL Genes (Not in any global winner) [{len(u_plus)}]:")
                for i in range(0, len(u_plus), 10): rprint(f"      {', '.join(u_plus[i:i+10])}")
            else:
                rprint(f"    * Methodology 100% contained within the Global Union.")
            if u_minus:
                label = f"[{len(u_minus)}]" if len(u_minus) < 100 else f"[{len(u_minus)}] (Truncated)"
                rprint(f"    * MISSING Union Genes {label}:")
                if len(u_minus) < 100:
                    for i in range(0, len(u_minus), 10): rprint(f"      {', '.join(u_minus[i:i+10])}")
                else:
                    rprint(f"      (Too many genes to list)")
        rprint("-" * 130)

    # Escribir el reporte al archivo
    with open(os.path.join(output_dir, "audit_full_report.txt"), "w", encoding="utf-8") as f:
        f.write(report.getvalue())
    
    rprint("\n" + "=" * 130)
    rprint(f"Audit completed. Results saved in: {output_dir}")
    rprint("=" * 130)

if __name__ == "__main__":
    run_methodology_audit()
