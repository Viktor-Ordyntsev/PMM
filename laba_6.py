import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Параметры задачи
a, b = 1, 1
Nx, Ny = 70, 70
dx, dy = a / (Nx - 1), b / (Ny - 1)

x = np.linspace(0, a, Nx)
y = np.linspace(0, b, Ny)
X, Y = np.meshgrid(x, y)

rho = 20 * (X**2 + Y**2)

u = np.zeros((Ny, Nx))

# Граничные условия
u[:, 0] = 3 * y 
u[:, -1] = 0     
u[0, :] = 0      
u[-1, :] = 3 * (1 - x**2) 

tolerance = 1e-6
max_iter = 10000
frames = []  

for it in range(max_iter):
    u_old = u.copy()
    
    for i in range(1, Ny-1):
        for j in range(1, Nx-1):
            u[i, j] = 0.25 * (
                u_old[i+1, j] + u_old[i-1, j] +
                u_old[i, j+1] + u_old[i, j-1] -
                dx**2 * rho[i, j]
            )
    
    if it % 5 == 0:  
        frames.append(u.copy())
    
    if np.max(np.abs(u - u_old)) < tolerance:
        print(f"Сошлось за {it} итераций.")
        frames.append(u.copy())
        break
else:
    print("Достигнуто максимальное число итераций.")

fig, ax = plt.subplots(figsize=(8, 6))
mesh = ax.pcolormesh(X, Y, frames[0], cmap='hot', shading='auto')
cbar = fig.colorbar(mesh)
cbar.set_label('u(x, y)')
ax.set_title("Решение уравнения Пуассона")
ax.set_xlabel("x")
ax.set_ylabel("y")
def update(frame):
    mesh.set_array(frame.ravel())
    return mesh,

ani = FuncAnimation(fig, update, frames=frames, interval=100, blit=True)
plt.show()