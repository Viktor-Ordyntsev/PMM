import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ПАРАМЕТРЫ ПРОГРАММЫ
L = 1.0                
N = 200                
h = L / (N - 1)        
tau = 0.000001           
T0 = 0                 
T_max = 1              
time_steps = int(T_max / tau)
sigma = 5.0            

x = np.linspace(0, L, N)
T = np.zeros(N)
T_new = np.zeros_like(T)


T[int(0.4 / h)] = T0

def tridiagonal_matrix_algorithm(A, B, C, D):
    n = len(D)
    P = np.zeros(n)
    Q = np.zeros(n)
    
    # прямой ход
    P[0] = -C[0] / B[0]
    Q[0] = D[0] / B[0]
    
    for i in range(1, n - 1):
        denom = B[i] + A[i] * P[i - 1]
        P[i] = -C[i] / denom
        Q[i] = (D[i] - A[i] * Q[i - 1]) / denom
    
    # обратный ход
    X = np.zeros(n)
    X[-1] = (D[-1] - A[-1] * Q[-2]) / (B[-1] + A[-1] * P[-2])
    
    for i in range(n - 2, -1, -1):
        X[i] = P[i] * X[i + 1] + Q[i]
    
    return X

# итерационный процесс неявной схемы
def implicit_scheme(T, N, h, tau, sigma):
    A = np.zeros(N - 1)
    B = np.zeros(N)
    C = np.zeros(N - 1)
    D = np.zeros(N)

    # заполнение коэффициентов для тридиагональной системы
    for j in range(1, N - 1):
        A[j - 1] = -sigma * tau / h**2  # нижняя диагональ
        B[j] = 1 + 2 * sigma * tau / h**2  # главная диагональ
        C[j] = -sigma * tau / h**2  # верхняя диагональ
        D[j] = T[j] + (1 - sigma) * tau / h**2 * (T[j + 1] - 2 * T[j] + T[j - 1])  # правые части

    
    D[0] = -1 #T[0]  # T(0, t) = 0
    D[-1] = 0 #T[-1]  # T(L, t) = 0
    B[0] = 1
    B[-1] = 1
    A[0] = 0
    C[-1] = 0
    
    # решение тридиагональной системы
    T_new = tridiagonal_matrix_algorithm(A, B, C, D)
    
    return T_new

# Сохранение температур для нескольких моментов времени
temperature_profiles = []
for time_step in range(0, time_steps, max(1, time_steps // 10)):
    for _ in range(time_steps // 10):
        T = implicit_scheme(T, N, h, tau, sigma)
    temperature_profiles.append(T.copy())

# Построение графиков
plt.figure(figsize=(10, 6))
time_steps_plot = np.arange(0, time_steps, max(1, time_steps // 10))  
for idx, time_step in enumerate(time_steps_plot):
    plt.plot(x, temperature_profiles[idx], label=f't = {time_step * tau:.2f} s')

plt.title('Решение уравнения теплопроводности - неявная схема')
plt.xlabel('x')
plt.ylabel('T(x,t)')
plt.legend()
plt.grid(True)
plt.show()