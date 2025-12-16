import matplotlib.pyplot as plt
import pandas as pd
import glob
import os

def plot_multi_k_comparison():
    # Configuración de los K que tienes disponibles
    k_targets = [15, 21, 31]
    
    # Configuración de colores para mantener consistencia visual
    chromosomes_meta = {
        '19': {'color': 'red', 'label': 'Chr 19'},
        '12': {'color': 'blue', 'label': 'Chr 12'},
        'Y': {'color': 'green', 'label': 'Chr Y'}
    }

    # Crear una figura con 3 subplots en una fila
    # sharey=True hace que todos compartan el eje Y para comparar magnitudes fácilmente
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    # Título general
    fig.suptitle('Evolución del Espectro de Frecuencia según longitud K', fontsize=16, y=1.05)

    # Iterar sobre cada K y asignarle un subplot
    for i, k in enumerate(k_targets):
        ax = axes[i]
        
        # Buscar archivos que coincidan con el patrón k específico dentro de la carpeta csv
        # Ajusta "csv/" si tus archivos están en otra ruta
        pattern = f"csv/*_k{k}_*.csv" 
        files = glob.glob(pattern)
        
        if not files:
            ax.text(0.5, 0.5, f"No se hallaron archivos\npara K={k}", 
                    ha='center', va='center', transform=ax.transAxes)
            continue

        for f in files:
            # Identificar de qué cromosoma es el archivo
            filename = os.path.basename(f)
            c_key = None
            
            # Lógica para detectar 12, 19 o MT en el nombre del archivo
            if "chromosome.19" in filename: c_key = '19'
            elif "chromosome.12" in filename: c_key = '12'
            elif "chromosome.Y" in filename: c_key = 'Y'
            
            if not c_key: continue # Saltar archivos desconocidos

            # Leer datos
            try:
                df = pd.read_csv(f)
            except Exception as e:
                print(f"Error leyendo {f}: {e}")
                continue
                
            df = df.sort_values("Frecuencia")

            # --- NORMALIZACIÓN (Clave para comparar) ---
            # Dividimos por el total para obtener Probabilidad
            total_kmers = df["Conteo"].sum()
            df["Probabilidad"] = df["Conteo"] / total_kmers

            # Plotear (Log-Log)
            ax.loglog(df["Frecuencia"], df["Probabilidad"], 
                      color=chromosomes_meta[c_key]['color'],
                      label=chromosomes_meta[c_key]['label'] if i == 0 else "", # Leyenda solo en el 1ro
                      linewidth=2, alpha=0.8)

        # Estética del Subplot
        ax.set_title(f"Longitud K = {k}", fontsize=14, fontweight='bold')
        ax.set_xlabel("Frecuencia (Veces que se repite)", fontsize=10)
        ax.grid(True, which="both", ls="--", alpha=0.3)
        
        # Solo poner etiqueta Y en el primer gráfico para no ensuciar
        if i == 0:
            ax.set_ylabel("Probabilidad (Densidad)", fontsize=12)

    # Crear una leyenda global arriba
    # Obtenemos los handles (líneas) y labels del primer ax
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.0), 
               ncol=3, fontsize=12, frameon=False)

    plt.tight_layout()
    plt.savefig("results/comparativa_multi_k.png", dpi=300, bbox_inches='tight')

def plot_scores():
    df = pd.read_csv('csv/resultados_scores.csv')

    # Convertimos a una lista de tuplas para mantener compatibilidad con tu lógica de ordenamiento
    # Forzamos 'cromosoma' a string para que '1' no sea entero todavía y funcione el sort_key
    data_list = list(zip(df['cromosoma'].astype(str), df['score']))

    # 2. Procesamiento de datos
    # MODIFICADO: Se eliminó el bucle de Regex porque el CSV ya tiene los nombres limpios ("1", "X", etc.)

    # Función para ordenar naturalmente (1, 2, ... 10, ... X, Y, MT)
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

    # 3. Configuración de Colores
    colors = []
    for n, s in zip(names, scores):
        if n == "MT":
            colors.append('#ff4d4d')  # Rojo para anomalía
        elif n == "19":
            colors.append('#2ecc71')  # Verde brillante para el más denso
        elif s > 40: # Nota: Con tus nuevos datos, solo el 19 supera 40, quizás quieras bajar este umbral a 35
            colors.append('#27ae60')  # Verde oscuro para otros muy densos
        elif n == "Y":
            colors.append('#f1c40f')  # Amarillo para el Y (bajo)
        else:
            colors.append('#3498db')  # Azul estándar

    # 4. Generar Gráfico
    plt.figure(figsize=(14, 7))
    bars = plt.bar(names, scores, color=colors, edgecolor='black', alpha=0.8)

    # Línea de base (0)
    plt.axhline(0, color='black', linewidth=1)

    # Etiquetas y Títulos
    plt.title('Score por Cromosoma', fontsize=16, fontweight='bold')
    plt.xlabel('Cromosoma', fontsize=12)
    plt.ylabel('Score', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    # Anotaciones en las barras clave
    for bar, name, score in zip(bars, names, scores):
        height = bar.get_height()
        # MODIFICADO: Agregué una pequeña validación para que el texto no se superponga si la barra es muy baja
        if name in ["19", "MT", "Y"]:
            y_pos = height + (1 if height > 0 else -3)
            # Ajuste específico para MT si es muy pequeño pero positivo
            if name == "MT" and 0 < height < 1: 
                y_pos = height + 0.5 
                
            plt.text(bar.get_x() + bar.get_width()/2., 
                    y_pos, 
                    f'{score:.1f}', 
                    ha='center', va='bottom' if height>0 else 'top', 
                    fontsize=10, fontweight='bold', color='black')

    # Leyenda personalizada
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#3498db', label='Estándar'),
        Patch(facecolor='#2ecc71', label='Alta Densidad (Cromosoma 19)'),
        Patch(facecolor='#ff4d4d', label='MT'),
        Patch(facecolor='#f1c40f', label='Baja Complejidad (Y)')
    ]
    plt.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig('results/resultados_score.png', dpi=300)

if __name__ == "__main__":
    plot_multi_k_comparison()
    plot_scores()
    print("Plots guardados en results/")