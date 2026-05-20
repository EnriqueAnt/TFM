import os
import pandas as pd
import glob

# Script para el análisis estratégico de ganadores basado en robustez estadística (SR) y rendimiento (Media)
def analyze_experiments():
    # Configuración de rutas
    # El script se ejecuta desde Project/info_extra/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    base_dir = os.path.join(project_root, "results_ml_classification")
    
    # Subcarpeta para resultados
    output_dir = os.path.join(script_dir, "winners_analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    report_txt_path = os.path.join(output_dir, "detailed_winners_report.txt")
    report_csv_path = os.path.join(output_dir, "summary_winners_table.csv")

    print("-" * 145)
    print(f"{'STRATEGIC WINNER ANALYSIS (SR VS MEAN)':^145}")
    print(f"Output directory: {output_dir}")
    print("-" * 145)

    if not os.path.exists(base_dir):
        print(f"Error: Results folder not found at {base_dir}")
        return

    # Buscar archivos de resumen
    pattern = os.path.join(base_dir, "*", "k_*", "best_analysis", "summary_signed_rank.csv")
    summary_files = glob.glob(pattern)

    if not summary_files:
        print("No results found in the results folder.")
        return

    # Sort files numerically by K (k_25, k_50, k_75, k_100...)
    def get_sort_key(path):
        p = path.split(os.sep)
        ds = p[-4]
        k_val = int(p[-3].split('_')[1])
        return (ds, k_val)
    
    summary_files.sort(key=get_sort_key)

    table_results = []
    
    with open(report_txt_path, "w", encoding="utf-8") as f_out:
        f_out.write("="*120 + "\n")
        f_out.write(f"{'REPORTE DETALLADO DE SELECCIÓN DE METODOLOGÍAS (ORDENADO POR K)':^120}\n")
        f_out.write("="*120 + "\n\n")

        for sr_path in summary_files:
            parts = sr_path.split(os.sep)
            dataset = parts[-4]
            k_val = parts[-3]
            mean_path = sr_path.replace("summary_signed_rank.csv", "summary_mean.csv")

            if not os.path.exists(mean_path): continue

            df_sr = pd.read_csv(sr_path, index_col=0)
            df_mean = pd.read_csv(mean_path, index_col=0)

            # Identificación de líderes de SR (Manejo de empates)
            max_sr = df_sr["Total"].max()
            ties_sr = df_sr[df_sr["Total"] == max_sr].index.tolist()
            sr_was_tied = len(ties_sr) > 1
            
            if sr_was_tied:
                # Excepción: Desempate por mejor posición en la Media
                winner_sr_name = df_mean.loc[ties_sr].sort_values(by="Total", ascending=False).index[0]
            else:
                winner_sr_name = ties_sr[0]

            # Líder de la Media
            winner_mean_name = df_mean.index[0]

            # Valores y puestos cruzados
            score_sr_in_sr = df_sr.loc[winner_sr_name, "Total"]
            score_sr_in_mean = df_mean.loc[winner_sr_name, "Total"]
            rank_sr_in_mean = list(df_mean.index).index(winner_sr_name) + 1

            score_mean_in_sr = df_sr.loc[winner_mean_name, "Total"]
            score_mean_in_mean = df_mean.loc[winner_mean_name, "Total"]
            rank_mean_in_sr = list(df_sr.index).index(winner_mean_name) + 1

            # Deltas y Fórmula
            d_sr = score_sr_in_sr - score_mean_in_sr
            d_m_over_sr = score_mean_in_mean - score_sr_in_mean
            d_m_rel = score_sr_in_mean - score_mean_in_mean # Negativo si el de media es mejor
            
            formula_val = (3 * d_sr) + d_m_rel
            do_match = (winner_sr_name == winner_mean_name) or (d_m_over_sr == 0)
            
            final_winner = winner_sr_name
            source = "SR"
            if not do_match and d_m_over_sr > 3 * d_sr:
                final_winner = winner_mean_name
                source = "Mean"

            # BA Analysis
            max_ba_sr_val = df_sr["BA"].max()
            best_ba_sr_methods = df_sr[df_sr["BA"] == max_ba_sr_val].index.tolist()
            winner_has_best_ba_sr = final_winner in best_ba_sr_methods

            max_ba_mean_val = df_mean["BA"].max()
            best_ba_mean_methods = df_mean[df_mean["BA"] == max_ba_mean_val].index.tolist()
            winner_has_best_ba_mean = final_winner in best_ba_mean_methods


            # Construcción del reporte en castellano para el TXT
            msg = f"DATASET: {dataset.upper()} | K: {k_val}\n" + "-"*60 + "\n"
            if do_match:
                msg += f"ESTADO: COINCIDENCIA (MATCH YES).\n"
                msg += f"Ganador: {final_winner}\n"
                msg += f"Score SR: {score_sr_in_sr} | Score Media: {score_sr_in_mean}\n"
                
                dist_sr_2nd = score_sr_in_sr - df_sr.iloc[1]["Total"] if len(df_sr) > 1 else 0
                dist_mean_2nd = score_sr_in_mean - df_mean.iloc[1]["Total"] if len(df_mean) > 1 else 0
                msg += f"Distancia al 2º en SR: {dist_sr_2nd} | Distancia al 2º en Media: {dist_mean_2nd}\n"
                
                if sr_was_tied:
                    msg += f"Aviso: Empate en liderazgo SR (Valor {max_sr}). Se eligió por mejor puesto en Media.\n"
                    msg += f"Candidatos empatados en SR:\n"
                    for t in ties_sr:
                        msg += f"  - {t}: Puesto Media #{list(df_mean.index).index(t)+1} ({df_mean.loc[t, 'Total']})\n"
            else:
                msg += f"ESTADO: CONFLICTO DE LIDERAZGO (MATCH NO).\n"
                msg += f"Líder SR: {winner_sr_name}\n"
                msg += f"  - Puesto en SR: #1 (Valor: {score_sr_in_sr})\n"
                msg += f"  - Puesto en Media: #{rank_sr_in_mean} (Valor: {score_sr_in_mean})\n"
                msg += f"Líder Media: {winner_mean_name}\n"
                msg += f"  - Puesto en Media: #1 (Valor: {score_mean_in_mean})\n"
                msg += f"  - Puesto en SR: #{rank_mean_in_sr} (Valor: {score_mean_in_sr})\n"
                
                msg += f"Cálculo de Decisión:\n"
                msg += f"  - Ventaja SR (dSR): {d_sr}\n"
                msg += f"  - Ventaja Media (dM): {d_m_over_sr}\n"
                msg += f"  - Valor Decisión (3*dSR + dM_rel): (3*{d_sr}) + ({d_m_rel}) = {formula_val}\n"
                
                if source == "Mean":
                    msg += f"RESULTADO: CAMBIO POR REGLA 3X. El líder de Media es el ganador final.\n"
                else:
                    msg += f"RESULTADO: PERSISTENCIA DEL LÍDER ESTADÍSTICO. No se triplica la ventaja.\n"
                
                # Metodologías intermedias
                idx_sr_in_mean = list(df_mean.index).index(winner_sr_name)
                intermediate = list(df_mean.index)[1:idx_sr_in_mean]
                if intermediate:
                    msg += f"Metodologías entre el Líder de Media y el Líder SR (Ranking Medias):\n"
                    for m in intermediate:
                        msg += f"  - {m}: Puesto Media #{list(df_mean.index).index(m)+1} ({df_mean.loc[m, 'Total']}), Puesto SR #{list(df_sr.index).index(m)+1} ({df_sr.loc[m, 'Total']})\n"

            msg += f"GANADOR FINAL: {final_winner}\n"
            
            msg += f"--- Análisis de Balanced Accuracy (BA) ---\n"
            if winner_has_best_ba_sr:
                msg += f"  - SR: El ganador final tiene el mayor BA en SR ({max_ba_sr_val}).\n"
            else:
                msg += f"  - SR: El ganador NO tiene el mayor BA en SR. Método(s) líder(es) en BA (SR): {', '.join(best_ba_sr_methods)} (Valor: {max_ba_sr_val})\n"
                
            if winner_has_best_ba_mean:
                msg += f"  - Media: El ganador final tiene el mayor BA en Media ({max_ba_mean_val}).\n"
            else:
                msg += f"  - Media: El ganador NO tiene el mayor BA en Media. Método(s) líder(es) en BA (Media): {', '.join(best_ba_mean_methods)} (Valor: {max_ba_mean_val})\n"
            msg += "\n"
            f_out.write(msg)

            # Registro de datos para el CSV
            table_results.append({
                "Dataset": dataset,
                "K": k_val,
                "Match": "YES" if do_match else "NO",
                "Source": source,
                "SR": df_sr.loc[final_winner, "Total"],
                "dSR": d_sr,
                "Mean": df_mean.loc[final_winner, "Total"],
                "dM": d_m_rel,
                "Decision": formula_val,
                "Winning Methodology": final_winner,
                "Best_BA_SR": "YES" if winner_has_best_ba_sr else "NO",
                "Best_BA_Mean": "YES" if winner_has_best_ba_mean else "NO"
            })

    df_final = pd.DataFrame(table_results)
    df_final.to_csv(report_csv_path, index=False, sep=";")

    print("\n" + "-"*145)
    print(f"{'SUMMARY TABLE OF FINAL WINNERS':^145}")
    print("-" * 145)
    
    console_header = f"{'Dataset':<10} | {'K':<7} | {'Match':<5} | {'Source':<6} | {'SR':<4} | {'dSR':<4} | {'Mean':<4} | {'dM':<6} | {'Decision':<8} | {'Best_BA_SR':<10} | {'Best_BA_Mean':<12} | {'Winning Methodology'}"
    print(console_header)
    print("-" * 145)
    for r in table_results:
        print(f"{r['Dataset']:<10} | {r['K']:<7} | {r['Match']:<5} | {r['Source']:<6} | "
              f"{r['SR']:<4.1f} | {r['dSR']:<4.1f} | {r['Mean']:<4.1f} | "
              f"{r['dM']:<6.1f} | {r['Decision']:<8.1f} | {r['Best_BA_SR']:<10} | {r['Best_BA_Mean']:<12} | {r['Winning Methodology']}")
    print("-" * 145)
    print(f"Analysis completed. Reports generated in: {output_dir}")

if __name__ == "__main__":
    analyze_experiments()