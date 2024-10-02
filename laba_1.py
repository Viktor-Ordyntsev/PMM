import numpy as np
import matplotlib.pyplot as plt

# Заданные параметры
L = 1.0  # длина стержня
T_max = 4.0  # конечное время
h = 0.02  # шаг по пространству
tau = 0.0001  # шаг по времени (для стабильности явной схемы)
x = np.arange(0, L + h, h)  # сетка по пространству
t = np.arange(0, T_max + tau, tau)  # сетка по времени
N_x = len(x)
N_t = len(t)

# Явная схема
alpha = tau / h**2
T_explicit = np.zeros((N_t, N_x))

# Начальные условия
T_explicit[0, :] = 0

# Функция для граничных условий на левой границе
def left_boundary(ti):
    return ti / (1 + ti)

# Явная схема
for n in range(0, N_t - 1):
    for i in range(1, N_x - 1):
        T_explicit[n + 1, i] = T_explicit[n, i] + alpha * (T_explicit[n, i - 1] - 2 * T_explicit[n, i] + T_explicit[n, i + 1])
    
    # Граничные условия
    T_explicit[n + 1, 0] = left_boundary(t[n + 1])  # T(0, t) = t / (1 + t)
    T_explicit[n + 1, -1] = 0  # T(1, t) = 0

# Построение графиков для нескольких моментов времени
plt.figure(figsize=(10, 6))
for time_step in [0, int(1/tau), int(2/tau), int(3/tau), int(4/tau)]:
    plt.plot(x, T_explicit[time_step, :], label=f't = {time_step * tau:.1f} s')

plt.title('Решение уравнения теплопроводности - явная схема')
plt.xlabel('x')
plt.ylabel('T(x,t)')
plt.legend()
plt.grid(True)
plt.show()