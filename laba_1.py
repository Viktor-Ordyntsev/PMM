import numpy as np
import matplotlib.pyplot as plt

L = 1       
T_max = 1  
alpha = 0.5     
N_x = 200        

h = L / (N_x - 1)         

tau = h**2 / (2 * alpha)
N_t = int(T_max / tau)

x = np.linspace(0, L, N_x)  
t = np.linspace(0, T_max, N_t)  


T_explicit = np.zeros((N_x, N_t))


T_explicit[:, 0] = 0  

for n in range(0, N_t - 1):
    for i in range(1, N_x - 1):
        T_explicit[i, n + 1] = ((tau * alpha) / (h * h)) * T_explicit[i+1, n] + (1 - 2 * ((tau * alpha) / (h * h))) * T_explicit[i, n] + ((tau * alpha) / (h * h)) * T_explicit[i-1, n]
    
    # Граничные условия
    T_explicit[0, n+1] = - 1  
    T_explicit[-1, n+1] = 0

# Построение графиков
plt.figure(figsize=(10, 6))
time_steps = np.arange(0, 0.55, 0.05)
for time_step in time_steps:
    index = int(time_step / tau)
    plt.plot(x, T_explicit[:, index], label=f't = {time_step:.2f} s')

plt.title('Решение уравнения теплопроводности - явная схема')
plt.xlabel('x')
plt.ylabel('T(x,t)')
plt.legend()
plt.grid(True)
plt.show()