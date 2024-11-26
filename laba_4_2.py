import numpy as np
import matplotlib.pyplot as plt

# ПАРАМЕТРЫ ПРОГРАММЫ
wind_speed = 1  # скорость ветра, м/с
D_min = 0.01  # минимальный коэффициент диффузии
D_max = 1.0  # максимальный коэффициент диффузии
C_min = 0.01  # минимальная скорость ветра для максимального коэффициента диффузии
D0 = 0.3  # коэффициент для расчета турбулентной диффузии
z_p = 50.0  # высота пограничного слоя

x_len = 20 # длина области по x, км
y_len = 20  # длина области по y, км
grid_size = 1000  # размер сетки

x = np.linspace(-x_len, x_len, grid_size)
y = np.linspace(-y_len, y_len, grid_size)
x, y = np.meshgrid(x, y)

def concentration(x, y, wind_speed, wind_direction, D_min, D_max, C_min, D0, z_p):
    Q = 5  # источник выброса (условная единица)
    
    # Поворот координат по направлению ветра
    x_rot = x * np.cos(wind_direction) + y * np.sin(wind_direction)
    y_rot = -x * np.sin(wind_direction) + y * np.cos(wind_direction)
    
    # Расчет коэффициентов и стандартных отклонений
    C_mod = np.maximum(np.sqrt(wind_speed**2), C_min)
    D = np.where(C_mod >= C_min, D_min + (D_max - D_min) * (C_mod / wind_speed), D_min)
    sigma_y = D * np.sqrt(np.maximum(x_rot / wind_speed, 1e-10))
    sigma_z = D0 * np.exp(-4.0 * x_rot / z_p) * (C_mod + 1)
    
    # Расчет концентрации
    C = (Q / (2 * np.pi * wind_speed * sigma_y * sigma_z)) * np.exp(-y_rot**2 / (2 * sigma_y**2))
    C[x_rot < 0] = 0  # Удаляем "хвосты" в обратном направлении

    return C

# Создание комбинированной модели
z_combined = np.zeros_like(x)

# Первая половина времени — ветер 45°
wind_direction_1 = np.radians(45)
z1 = concentration(x, y, wind_speed, wind_direction_1, D_min, D_max, C_min, D0, z_p)
z_combined += z1

# Вторая половина времени — ветер 315°
wind_direction_2 = np.radians(315)
z2 = concentration(x, y, wind_speed, wind_direction_2, D_min, D_max, C_min, D0, z_p)
z_combined += z2

# Итоговая визуализация
plt.figure(figsize=(8, 8))
cs = plt.contour(x, y, z_combined, levels=np.logspace(-3, 0, 10), colors='black')
plt.clabel(cs, inline=1, fontsize=10)
plt.title('Диффузионный перенос примесей (с учетом смены ветра)')
plt.xlabel('X (км)')
plt.ylabel('Y (км)')
plt.xlim(-x_len, y_len)
plt.ylim(-x_len, y_len)
plt.grid(True)
plt.show()