import numpy as np
import matplotlib.pyplot as plt

# Параметры задачи
h = 0.02  # Шаг по пространству
tau = 0.0001  # Шаг по времени
L = 1  # Длина стержня
T_max = 2 # Время моделирования

# Число шагов по пространству и времени
N = int(L / h)
M = int(T_max / tau)

# Коэффициент схемы
r = tau / h**2

# Инициализация сетки
T = np.zeros((M+1, N+1))

# Начальные условия T(x, 0) = 0 уже заданы нулями

# Итерации по времени
for n in range(M):
    # Граничное условие в точке x=1: T(1, t) = 0
    T[n+1, N] = 0
    
    # Граничное условие T_x(0, t) = -1 -> T_0 = T_1 + h
    T[n+1, 0] = T[n+1, 1] + h
    
    # Основное уравнение для внутренних точек
    for i in range(1, N):
        T[n+1, i] = T[n, i] + r * (T[n, i+1] - 2*T[n, i] + T[n, i-1])

# Построение графика температурного распределения в разные моменты времени
x = np.linspace(0, L, N+1)

plt.figure(figsize=(10,6))

time_steps = [0, M//10, M//5, M//2, M]  # Моменты времени для графика

for n in time_steps:
    plt.plot(x, T[n, :], label=f't={n*tau:.2f} s')

plt.xlabel('x')
plt.ylabel('T(x, t)')
plt.legend()
plt.title('Распределение температуры вдоль стержня на различных временных шагах')
plt.grid(True)
plt.show()