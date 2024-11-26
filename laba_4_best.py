import numpy as np
import matplotlib.pyplot as plt

# ПАРАМЕТРЫ ПРОГРАММЫ
wind_speed = 0.5  # скорость ветра, м/с
D_min = 0.1  # минимальный коэффициент диффузии
D_max = 1.0  # максимальный коэффициент диффузии
C_min = 0.01  # минимальная скорость ветра для максимального коэффициента диффузии
D0 = 0.3  # коэффициент для расчета турбулентной диффузии
z_p = 50.0  # высота пограничного слоя

x_len = 80  # длина области по x, км
y_len = 80  # длина области по y, км
grid_size = 1000  # размер сетки

x = np.linspace(-x_len, x_len, grid_size)
y = np.linspace(-y_len, y_len, grid_size)
x, y = np.meshgrid(x, y)

def concentration(x, y, wind_speed, wind_direction, D_min, D_max, C_min, D0, z_p):
    Q = 1  # источник выброса (условная единица)
    
    x_rot = x * np.cos(wind_direction) + y * np.sin(wind_direction)
    y_rot = -x * np.sin(wind_direction) + y * np.cos(wind_direction)
    
    C_mod = np.maximum(np.sqrt(wind_speed**2), C_min)
    D = np.where(C_mod >= C_min, D_min + (D_max - D_min) * (C_mod / wind_speed), D_min)
    
    sigma_y = D * np.sqrt(np.maximum(x_rot / wind_speed, 1e-10))
    sigma_z = D0 * np.exp(-4.0 * x_rot / z_p) * (C_mod + 1)
    
    C = (Q / (2 * np.pi * wind_speed * sigma_y * sigma_z)) * np.exp(-y_rot**2 / (2 * sigma_y**2))
    
    C = (Q / (2 * np.pi * wind_speed * sigma_y * sigma_z)) * np.exp(-y_rot**2 / (2 * sigma_y**2))
    C[x_rot < 0] = 0

    # Применение граничных условий
    C[x <= -x_len] = 0
    C[x >= x_len] = 0
    C[y <= -y_len] = 0
    C[y >= y_len] = 0
    
    return C

# Основной расчет
angles = [np.pi / 4, 7 * np.pi / 4]  # Углы: 45 и 315 градусов
z = np.zeros_like(x)

for i, angle in enumerate(angles):
    z += concentration(x, y, wind_speed, angle, D_min, D_max, C_min, D0, z_p) * (0.5 if i == 0 else 0.5)

# Построение графика
plt.figure(figsize=(8, 8))
cs = plt.contour(x, y, z, levels=np.logspace(-3, 0, 10), colors='black')
plt.clabel(cs, inline=1, fontsize=10)
plt.title('Диффузионный перенос примесей')
plt.xlabel('X (км)')
plt.ylabel('Y (км)')
plt.xlim(-x_len, x_len)
plt.ylim(-y_len, y_len)
plt.grid(True)
plt.show()