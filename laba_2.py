import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Параметры задачи
L = 1.0       
T_max = 1  
alpha = 0.5     
N = 200         
       
h = L / (N - 1)         
tau = h**2 / (2 * alpha)

M = int(T_max / tau)


x = np.linspace(0, L, N)
t = np.linspace(0, T_max, M)


T = np.zeros((M, N))


T[0, :] = 0


for n in range(M - 1):
    
    A = np.zeros(N)
    B = np.zeros(N)
    
    
    A[1] = 0
    B[1] = -1 
    
    for i in range(1, N - 1):
        A[i + 1] = (tau*alpha) / (h*h + (tau*alpha) * (2 - A[i]))
        B[i + 1] = (h*h * T[n, i] + (tau*alpha)* B[i] ) / (h*h + (tau*alpha) * (2 - A[i]))
        
    
    
    T_next = np.zeros(N)
    T_next[0] = -1

    
    for i in range(N - 2, 0, -1):
        T_next[i] = A[i+1] * T_next[i + 1] + B[i+1]
    
    
    T_next = np.maximum(T_next, -1) 
    T[n + 1, :] = T_next


# fig = plt.figure(figsize=(12, 6))


# ax1 = fig.add_subplot(121, projection='3d')
# X, Y = np.meshgrid(x,t)
# ax1.plot_surface(X, Y, T, cmap='hot')
# ax1.set_xlabel("Положение x вдоль стержня")
# ax1.set_ylabel("Время, t")
# ax1.set_zlabel("Температура")
# ax1.set_title("Распределение температуры (3D)")

# ax2 = fig.add_subplot(122)
# cax = ax2.imshow(T, extent=[0, L, 0, T_max], aspect='auto', origin='lower', cmap='hot')
# ax2.set_xlabel("Положение x вдоль стержня")
# ax2.set_ylabel("Время, t")
# ax2.set_title("Распределение температуры (2D)")
# fig.colorbar(cax, ax=ax2, label="Температура")

# def update(frame):
#     cax.set_array(T[:frame, :])
#     return cax,

# ani = FuncAnimation(fig, update, frames=M, interval=50, blit=True)

# plt.tight_layout()
# plt.show()

plt.figure(figsize=(10, 6))
time_steps = np.arange(0, 0.55, 0.05)
for time_step in time_steps:
    index = int(time_step / tau)
    plt.plot(x, T[index, :], label=f't = {time_step:.2f} s')

plt.title('Решение уравнения теплопроводности - НЕявная схема')
plt.xlabel('x')
plt.ylabel('T(x,t)')
plt.legend()
plt.grid(True)
plt.show()