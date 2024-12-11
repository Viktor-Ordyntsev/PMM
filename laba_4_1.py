import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Параметры модели
Lx, Ly = 20, 20 
Nx, Ny = 100, 100   
dx, dy = Lx / Nx, Ly / Ny  
dt = 0.01  
T = 5  
D = 0.5  # Коэффициент диффузии
alpha = 0.01  # Параметр вымывания
Q = 10  # Источник примеси

rho = np.zeros((Nx, Ny))
sources = [(10, 10)]  

# Инициализация начального состояния с учетом источников
for sx, sy in sources:
    ix, iy = int(sx / dx), int(sy / dy)  # Индексы в массиве
    rho[ix, iy] = Q  # Установить максимальную концентрацию в источнике


def boundary_conditions(rho):
    rho[0, :] = rho[1, :]
    rho[-1, :] = rho[-2, :]
    rho[:, 0] = rho[:, 1]
    rho[:, -1] = rho[:, -2]
    return rho


def update_rho(rho, dt, D, Cx, Cy, alpha, Q, sources):
    rho_new = rho.copy()
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            # Диффузия
            diffusion = D * (
                (rho[i+1, j] - 2*rho[i, j] + rho[i-1, j]) / dx**2 +
                (rho[i, j+1] - 2*rho[i, j] + rho[i, j-1]) / dy**2
            )
            # Адвекция
            advection = (
                Cx * (rho[i, j] - rho[i-1, j]) / dx +
                Cy * (rho[i, j] - rho[i, j-1]) / dy
            )
            # Источник и вымывание
            source_sink = -alpha * rho[i, j]
            for sx, sy in sources:
                if (i * dx - sx)**2 + (j * dy - sy)**2 < 0.1:
                    source_sink += Q
            # Обновление значения
            rho_new[i, j] = rho[i, j] + dt * (diffusion - advection + source_sink)
    return rho_new



fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(rho, extent=[0, Lx, 0, Ly], origin='lower', cmap='inferno', vmin=0, vmax=Q/4)
fig.colorbar(im, label='Плотность примеси')
ax.set_xlabel('X, км')
ax.set_ylabel('Y, км')
ax.set_title('Распределение примеси в атмосфере')



def animate(frame):
    global rho
    if frame < frames // 2:
        Cx_local, Cy_local = -2, 2
    else:
        Cx_local, Cy_local = -2, 2
    
    rho = update_rho(rho, dt, D, Cx_local, Cy_local, alpha, Q, sources)
    rho = boundary_conditions(rho)
    im.set_array(rho)
    return [im]

frames = int(T / dt)
ani = FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)

plt.show()
