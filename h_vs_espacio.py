import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Cargar los resultados del CSV generado por Fortran
df = pd.read_csv('resultados_mac_cormack_trapecio.csv', sep=';')

# 2. Tiempos objetivo del libro (Ejemplo 14-1)
tiempos_objetivo = [0, 500, 1000, 1500, 2000]

sns.set_theme(style="whitegrid")
plt.figure(figsize=(9, 6))

# 3. Graficar filtrando por el 'paso' único más cercano
for t_obj in tiempos_objetivo:
    # Obtener el índice del registro cuyo tiempo esté más cerca del objetivo
    idx_cercano = (df['tiempo'] - t_obj).abs().idxmin()
    
    # Obtener el número de paso numérico único (entero)
    paso_cercano = df.loc[idx_cercano, 'paso']
    t_real = df.loc[idx_cercano, 'tiempo']
    
    # Filtrar estrictamente por el número de PASO (evita duplicados de tiempo)
    subset = df[df['paso'] == paso_cercano].drop_duplicates(subset=['x']).sort_values('x')
    
    # Graficar h(x)
    plt.plot(subset['x'] / 1000.0, subset['h'], linewidth=2, label=f't = {int(round(t_real))} s')

# 4. Configuración de ejes y títulos
plt.title('Perfil espacial de la lámina de agua h(x)', fontsize=13, fontweight='bold')
plt.xlabel('Distancia (km)', fontsize=12)
plt.ylabel('Profundidad del flujo h (m)', fontsize=12)
plt.xlim(0, 5)
plt.ylim(5.3, 7.1)
plt.legend(title='Instante de Tiempo', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()