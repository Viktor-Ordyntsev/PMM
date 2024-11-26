import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Параметры задачи
L = 1.0        
T_max = 4.0    
h = 0.02       
tau = 0.005    
N = int(L / h) + 1  
M = int(T_max / tau) + 1  


x = np.linspace(0, L, N)
t = np.linspace(0, T_max, M)


u = np.zeros((M, N))

def c(x):
    return 1.0 if x <= 0.5 else 0.5


for i in range(N):
    if x[i] <= 1/3:
        u[0, i] = np.sin(3 * np.pi * x[i])**2
        u[1, i] = u[0, i] - tau * 3 * np.pi * np.sin(3 * np.pi * x[i])
    else:
        u[0, i] = 0
        u[1, i] = 0


for n in range(1, M - 1):
    for i in range(1, N - 1):
        
        c_i_plus_half = (c(x[i]) + c(x[i + 1])) / 2
        c_i_minus_half = (c(x[i]) + c(x[i - 1])) / 2
        
        
        u[n + 1, i] = (2 * u[n, i] - u[n - 1, i]
                       + (tau ** 2 / h ** 2) * (
                           c_i_plus_half * (u[n, i + 1] - u[n, i])
                           - c_i_minus_half * (u[n, i] - u[n, i - 1])))

   
    u[n + 1, 0] = 0
    u[n + 1, -1] = 0


fig, ax = plt.subplots()
ax.set_xlim(0, L)
ax.set_ylim(-1, 1)
line, = ax.plot(x, u[0, :], color='red')
ax.set_title("Волновое уравнение")
ax.set_xlabel("Позиция Х")
ax.set_ylabel("u(x, t)")


def animate(n):
    line.set_ydata(u[n, :])
    return line,


ani = FuncAnimation(fig, animate, frames=M, interval=50, blit=True)
plt.show()