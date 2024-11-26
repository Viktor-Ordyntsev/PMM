import numpy as np
import matplotlib.pyplot as plt

# Параметры
Lx, Ly = 80, 80
Nx, Ny = 1000, 1000
dx, dy = Lx / Nx, Ly / Ny
D = 0.2  # Коэффициент диффузии
alpha = 0.01  # Параметр вымывания
Q = 20  # Источник примеси
dt = 0.01  
T = 5  # Время моделирования

# Начальная концентрация (вся сетка пуста)
rho = np.zeros((Nx, Ny))

# Источник загрязнения (в центре)
sources = [(400, 400)]

# Граничные условия (нет потоков на границах)
def boundary_conditions(rho):
    rho[0, :] = rho[1, :]
    rho[-1, :] = rho[-2, :]
    rho[:, 0] = rho[:, 1]
    rho[:, -1] = rho[:, -2]
    return rho

# Вычисление компонент скорости ветра из угла
def wind_components(C, theta):
    Cx = C * np.cos(theta)  # Компонента по X
    Cy = C * np.sin(theta)  # Компонента по Y
    return Cx, Cy

# Обновление концентрации примеси
def update_rho(rho, dt, D, Cx, Cy, alpha, Q, sources, dx, dy):
    rho_new = np.copy(rho)
    
    # Вычисляем диффузию
    diffusion = D * (
        (np.roll(rho, -1, axis=0) - 2 * rho + np.roll(rho, 1, axis=0)) / dx**2 +
        (np.roll(rho, -1, axis=1) - 2 * rho + np.roll(rho, 1, axis=1)) / dy**2
    )
    
    # Вычисляем адвекцию
    advection_x = -Cx * (np.roll(rho, -1, axis=0) - rho) / dx
    advection_y = -Cy * (np.roll(rho, -1, axis=1) - rho) / dy
    advection = advection_x + advection_y
    
    # Источники и поглощение
    source_sink = -alpha * rho
    for sx, sy in sources:
        mask = (np.arange(Nx)[:, None] - sx / dx)**2 + (np.arange(Ny)[None, :] - sy / dy)**2 < (0.5 / dx)**2
        source_sink[mask] += Q
    
    rho_new += dt * (diffusion + advection + source_sink)
    return rho_new

# Параметры ветра
C = 5  # Скорость ветра (м/с)

# Моделируем процесс для времени от 0 до T
for t in range(int(T / dt)):
    # Меняем угол ветра со временем
    if t * dt < T / 2:
        theta = (45 * np.pi / 180)
    else:
         theta = (45 * np.pi / 180) #(135 * np.pi / 180)
    
    # Вычисляем компоненты скорости ветра
    Cx, Cy = wind_components(C, theta)
    
    # Обновляем концентрацию примеси
    rho = update_rho(rho, dt, D, Cx, Cy, alpha, Q, sources, dx, dy)
    
    # Применяем граничные условия
    rho = boundary_conditions(rho)

plt.figure(figsize=(8, 8))
cs = plt.contour(dx, dy, rho, levels=np.logspace(-3, 0, 10), colors='black')
plt.clabel(cs, inline=1, fontsize=10)
plt.title('Диффузионный перенос примесей')
plt.xlabel('X (км)')
plt.ylabel('Y (км)')
plt.xlim(-20, 20)
plt.ylim(-20, 20)
plt.grid(True)
plt.show()