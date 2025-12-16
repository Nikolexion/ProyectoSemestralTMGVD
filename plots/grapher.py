import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
from matplotlib.patches import Patch

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

CSV_DIR = os.path.join(SCRIPT_DIR, "csv")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")

def plot_multi_k_comparison():
    k_targets = [15, 21, 31]
    
    chromosomes_meta = {
        '19': {'color': 'red', 'label': 'Chr 19'},
        '12': {'color': 'blue', 'label': 'Chr 12'},
        'Y': {'color': 'green', 'label': 'Chr Y'}
    }

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    fig.suptitle('Evolución del Espectro de Frecuencia según longitud K', fontsize=16, y=1.05)

    for i, k in enumerate(k_targets):
        ax = axes[i]
        
        pattern = os.path.join(CSV_DIR, f"*_k{k}_*.csv")
        files = glob.glob(pattern)
        
        if not files:
            ax.text(0.5, 0.5, f"No se hallaron archivos\npara K={k}", 
                    ha='center', va='center', transform=ax.transAxes)
            continue

        for f in files:
            filename = os.path.basename(f)
            c_key = None
            
            if "chromosome.19" in filename: c_key = '19'
            elif "chromosome.12" in filename: c_key = '12'
            elif "chromosome.Y" in filename: c_key = 'Y'
            
            if not c_key: continue 

            try:
                df = pd.read_csv(f)
            except Exception as e:
                print(f"Error leyendo {f}: {e}")
                continue
                
            df = df.sort_values("Frecuencia")

            total_kmers = df["Conteo"].sum()
            df["Probabilidad"] = df["Conteo"] / total_kmers

            ax.loglog(df["Frecuencia"], df["Probabilidad"], 
                      color=chromosomes_meta[c_key]['color'],
                      label=chromosomes_meta[c_key]['label'] if i == 0 else "", #
                      linewidth=2, alpha=0.8)

        ax.set_title(f"Longitud K = {k}", fontsize=14, fontweight='bold')
        ax.set_xlabel("Frecuencia (Veces que se repite)", fontsize=10)
        ax.grid(True, which="both", ls="--", alpha=0.3)
        
        if i == 0:
            ax.set_ylabel("Probabilidad (Densidad)", fontsize=12)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.0), 
               ncol=3, fontsize=12, frameon=False)

    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "comparativa_multi_k.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

def plot_scores():

    csv_path = os.path.join(CSV_DIR, 'resultados_scores.csv')
    
    if not os.path.exists(csv_path):
        print(f"No se encuentra el archivo de scores en: {csv_path}")
        return

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"No se pudo leer el CSV: {e}")
        return

    data_list = list(zip(df['Archivo'].astype(str), df['Score']))

    def sort_key(item):
        key = item[0]
        if key.isdigit():
            return int(key)
        if key == "X": return 23
        if key == "Y": return 24
        if key == "MT": return 25
        return 26

    data_list.sort(key=sort_key)

    names = [x[0] for x in data_list]
    scores = [x[1] for x in data_list]

    colors = ['#3498db']

    plt.figure(figsize=(14, 7))
    bars = plt.bar(names, scores, color=colors, edgecolor='black', alpha=0.8)

    plt.axhline(0, color='black', linewidth=1)

    plt.title('Score por Archivo', fontsize=16, fontweight='bold')
    plt.xlabel('Archivo', fontsize=12)
    plt.ylabel('Score', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    for bar, name, score in zip(bars, names, scores):
        height = bar.get_height()

        if name in ["19", "MT", "Y"]:
            y_pos = height + (1 if height > 0 else -3)
            if name == "MT" and 0 < height < 1: 
                y_pos = height + 0.5 
                
            plt.text(bar.get_x() + bar.get_width()/2., 
                    y_pos, 
                    f'{score:.1f}', 
                    ha='center', va='bottom' if height>0 else 'top', 
                    fontsize=10, fontweight='bold', color='black')

    output_path = os.path.join(RESULTS_DIR, 'resultados_score.png')
    plt.savefig(output_path, dpi=300)

if __name__ == "__main__":
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR, exist_ok=True)
    
    plot_multi_k_comparison()
    plot_scores()
    print("Plots guardados en plots/results/")