import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 1. Cargar los resultados generados por solver_mac_cormack.90
df = pd.read_csv('resultados_mac_cormack_trapecio.csv', sep=';')

# Estilo visual de la gráfica
sns.set_theme(style="whitegrid")
plt.figure(figsize=(9, 6))

# 2. Selección explícita de nodos según el numeral 3 del libro
# Se requiere graficar la variación de h(t) a las distancias de 1.5, 2.5, 3.0 y 5.0 km
distancias_km = [1.5, 2.5, 3.0, 5.0]
distancias_m = [d * 1000 for d in distancias_km] # Convertimos a metros para buscar en la columna 'x'

nodos_seleccionados = []
for dist_km, dist_m in zip(distancias_km, distancias_m):
    # Calcular la diferencia absoluta para encontrar el nodo más cercano a la distancia objetivo
    diferencias = np.abs(df['x'] - dist_m)
    idx_min = diferencias.idxmin()
    nodo_cercano = df.loc[idx_min, 'nodo']
    nodos_seleccionados.append((nodo_cercano, dist_km))

colores = plt.cm.tab10(np.linspace(0, 0.6, len(nodos_seleccionados)))

# 3. Graficar h(t) para los nodos seleccionados
for (nodo, dist_km), color in zip(nodos_seleccionados, colores):
    subset = df[df['nodo'] == nodo]
    pos_x = subset['x'].iloc[0]
    
    # Se etiqueta indicando la distancia teórica en km y la ubicación x real del nodo evaluado
    plt.plot(subset['tiempo'], subset['h'], 
             linewidth=2, 
             label=f'x = {dist_km} km (Nodo {int(nodo)}, real: {pos_x:.1f} m)', 
             color=color)

# Títulos y etiquetas de los ejes
plt.title('Evolución temporal de la lámina de Agua h(t)', fontsize=14, fontweight='bold', pad=12)
plt.xlabel('Tiempo t (s)', fontsize=12, labelpad=8)
plt.ylabel('Profundidad del flujo h (m)', fontsize=12, labelpad=8)
plt.xlim(0, 2000)
plt.xticks(np.arange(0, 2001, 200))
plt.legend(title='Ubicación en el Canal', title_fontsize='11', fontsize=10, loc='best')
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()