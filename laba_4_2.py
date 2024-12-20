import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ПАРАМЕТРЫ ПРОГРАММЫ
wind_speed = 0.5  # скорость ветра, м/с
D_min = 0.1  # минимальный коэффициент диффузии
D_max = 1.0  # максимальный коэффициент диффузии
C_min = 0.01  # минимальная скорость ветра для максимального коэффициента диффузии
D0 = 0.3  # коэффициент для расчета турбулентной диффузии
z_p = 50.0  # высота пограничного слоя

x_len = 80  # длина области по x, км
y_len = 80  # длина области по y, км
grid_size = 300  # уменьшенный размер сетки для ускорения анимации

x = np.linspace(-x_len, x_len, grid_size)
y = np.linspace(-y_len, y_len, grid_size)
x, y = np.meshgrid(x, y)

def concentration(x, y, wind_speed, wind_direction, D_min, D_max, C_min, D0, z_p):
    Q = 1  # интенсивность выброса
    
    x_rot = x * np.cos(wind_direction) + y * np.sin(wind_direction)
    y_rot = -x * np.sin(wind_direction) + y * np.cos(wind_direction)
    
    C_mod = np.maximum(np.sqrt(wind_speed**2), C_min)
    D = np.where(C_mod >= C_min, D_min + (D_max - D_min) * (C_mod / wind_speed), D_min)
    
    sigma_y = D * np.sqrt(np.maximum(x_rot / wind_speed, 1e-10))
    sigma_z = D0 * np.exp(-4.0 * x_rot / z_p) * (C_mod + 1)
    
    C = (Q / (2 * np.pi * wind_speed * sigma_y * sigma_z)) * np.exp(-y_rot**2 / (2 * sigma_y**2))
    C[x_rot < 0] = 0

    # Применение граничных условий
    C[x <= -x_len] = 0
    C[x >= x_len] = 0
    C[y <= -y_len] = 0
    C[y >= y_len] = 0
    
    return C

# Настройка анимации
fig, ax = plt.subplots(figsize=(8, 8))
levels = np.logspace(-3, 0, 10)
contour = ax.contourf(x, y, np.zeros_like(x), levels=levels, cmap='viridis')
ax.set_title('Анимация перемешивания примесей')
ax.set_xlabel('X (км)')
ax.set_ylabel('Y (км)')
ax.set_xlim(-x_len, x_len)
ax.set_ylim(-y_len, y_len)
ax.grid(True)

angles = [np.pi / 2, 7 * np.pi / 4]  # Углы: 45 и 315 градусов

def update(frame):
    wind_direction = angles[frame % len(angles)]
    z = concentration(x, y, wind_speed, wind_direction, D_min, D_max, C_min, D0, z_p)
    for c in contour.collections:
        c.remove()  # Удаление предыдущих контуров
    return ax.contourf(x, y, z, levels=levels, cmap='viridis')

# Анимация
ani = FuncAnimation(fig, update, frames=20, interval=200, repeat=True)

plt.show()