import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import matplotlib.animation as animation

# 1. Cargar datos del CSV generado por Fortran
df = pd.read_csv('resultados_mac_cormack_trapecio.csv', sep=';')

# Filtrar y ordenar datos asegurando pasos únicos
df = df.sort_values(by=['paso', 'x']).drop_duplicates(subset=['paso', 'x'])

pasos_unicos = df['paso'].unique()
tiempos_dict = df.groupby('paso')['tiempo'].first().to_dict()

# 2. Configurar la figura y los ejes
fig, ax = plt.subplots(figsize=(10, 6))
plt.subplots_adjust(bottom=0.25)  # Espacio reservado para los controles abajo

x_km = df['x'] / 1000.0
h_min, h_max = df['h'].min(), df['h'].max()
padding = (h_max - h_min) * 0.1

ax.set_xlim(0, x_km.max())
ax.set_ylim(h_min - padding, h_max + padding)

ax.set_title('Evolución Temporal de la Lámina de Agua h(x, t) [MacCormack]', fontsize=13, fontweight='bold')
ax.set_xlabel('Distancia x (km)', fontsize=11)
ax.set_ylabel('Profundidad / Lámina de Agua h (m)', fontsize=11)
ax.grid(True, linestyle='--', alpha=0.7)

# Estado inicial (paso 0)
paso_inicial = pasos_unicos[0]
data_init = df[df['paso'] == paso_inicial]
line, = ax.plot(data_init['x'] / 1000.0, data_init['h'], color='#1f77b4', lw=2.5, label='Lámina de agua h(x)')

time_text = ax.text(0.02, 0.92, f'Tiempo: {tiempos_dict[paso_inicial]:.1f} s', 
                    transform=ax.transAxes, fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))

ax.legend(loc='upper right')

# 3. Crear controles interactivos (Slider y Botón)
ax_slider = plt.axes([0.15, 0.1, 0.55, 0.03], facecolor='lightgoldenrodyellow')
slider_paso = Slider(
    ax=ax_slider,
    label='Paso / Tiempo ',
    valmin=0,
    valmax=len(pasos_unicos) - 1,
    valinit=0,
    valfmt='%d'
)

ax_button = plt.axes([0.75, 0.095, 0.12, 0.04])
button = Button(ax_button, 'Pause', color='lightgray', hovercolor='0.9')

# Estado de animación
anim_running = True

def update(val):
    idx = int(slider_paso.val)
    paso = pasos_unicos[idx]
    sub = df[df['paso'] == paso]
    line.set_ydata(sub['h'])
    time_text.set_text(f'Tiempo: {tiempos_dict[paso]:.1f} s')
    fig.canvas.draw_idle()

slider_paso.on_changed(update)

def animate(frame):
    if anim_running:
        current_val = slider_paso.val
        next_val = (current_val + 1) % len(pasos_unicos)
        slider_paso.set_val(next_val)
    return line, time_text

def toggle_anim(event):
    global anim_running
    anim_running = not anim_running
    button.label.set_text('Play' if not anim_running else 'Pause')

button.on_clicked(toggle_anim)

# Iniciar animación (intervalo de 40 ms entre fotogramas)
anim = animation.FuncAnimation(fig, animate, interval=40, blit=False)

plt.show()