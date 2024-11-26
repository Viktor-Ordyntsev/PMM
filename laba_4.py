import numpy as np
import matplotlib.pyplot as plt


Lx, Ly = 50, 50 
Nx, Ny = 100, 100   
dx, dy = Lx / Nx, Ly / Ny  
dt = 0.01  
T = 5  
D = 0.2  # Коэффициент диффузии
Cx, Cy = 0, 1  # Скорость ветра (км/ч)
alpha = 0.01  # Параметр вымывания
Q = 10  # Источник примеси


rho = np.zeros((Nx, Ny))


sources = [(25, 25)]  # Координаты источника загрязнения

# Граничные (нет потоков на границах)
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
          
            diffusion = D * (
                (rho[i+1, j] - 2*rho[i, j] + rho[i-1, j]) / dx**2 +
                (rho[i, j+1] - 2*rho[i, j] + rho[i, j-1]) / dy**2
            )
            
            advection = (
                Cx * (rho[i, j] - rho[i-1, j]) / dx +
                Cy * (rho[i, j] - rho[i, j-1]) / dy
            )
            
            source_sink = -alpha * rho[i, j]
            for sx, sy in sources:
                if (i * dx - sx)**2 + (j * dy - sy)**2 < 0.1:
                    source_sink += Q
            
            rho_new[i, j] = rho[i, j] + dt * (diffusion - advection + source_sink)
    return rho_new


for t in range(int(T / dt)):
    rho = update_rho(rho, dt, D, Cx, Cy, alpha, Q, sources)
    rho = boundary_conditions(rho)





plt.figure(figsize=(8, 8))
plt.imshow(rho, extent=[0, Lx, 0, Ly], origin='lower', cmap='viridis')
plt.colorbar(label='Плотность примеси')
plt.xlabel('X, км')
plt.ylabel('Y, км')
plt.title('Распределение примеси в атмосфере')
plt.show()