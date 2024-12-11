import numpy as np
import matplotlib.pyplot as plt

# ПАРАМЕТРЫ ПРОГРАММЫ  
kappa_0 = 0.5            # начальная теплопроводность
u = 1                  # скорость волны
beta = 2            # степень нелинейности
xi_0 = 0.5

L = 1.0                
T_max = 1.0           
Nx = 200               
Nt = 10000              
dx = L / Nx            
global dt
dt = T_max / Nt        

x = np.linspace(0, L, Nx)
T = np.zeros(Nx)


T[int(xi_0 * Nx)] = 0


def T_xi(xi, t):
    #return (u * beta / kappa_0) ** (1 / beta) * (xi_0 - xi * u * t) ** (1 / beta)
    return kappa_0 * (t ** beta)


T_explicit = np.zeros((Nx, Nt))


for frame in range(1, Nt):
    T_new = np.zeros_like(T)
    
    for i in range(1, Nx - 1):
        T_new[i] = T[i] + T_xi(i / Nx, frame * dt) * (T[i + 1] - 2 * T[i] + T[i - 1])

    
    T_new[0] = -1   # при x = 0
    T_new[-1] = 0   # при x = 1

    T[:] = T_new[:]
    T_explicit[:, frame] = T[:]  

# Построение графиков
plt.figure(figsize=(10, 6))
time_steps = np.arange(0, Nt, int(Nt / 10))  
for time_step in time_steps:
    plt.plot(x, T_explicit[:, time_step], label=f't = {time_step * dt:.2f} s')

plt.title('Решение уравнения теплопроводности - явная схема')
plt.xlabel('x')
plt.ylabel('T(x,t)')
plt.legend()
plt.grid(True)
plt.show()