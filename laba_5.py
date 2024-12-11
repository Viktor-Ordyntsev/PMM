import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ПАРАМЕТРЫ ПРОГРАММЫ
kappa_0 = 5           # начальная теплопроводность
u = 2                  # скорость волны
beta = 2.0            # степень нелинейности
xi_0 = 0.2             # начальное положение волны

L = 1.0                
T_max = 1.0            
Nx = 100               
Nt = 10000             
dx = L / Nx            
global dt
dt = T_max / Nt        

# Инициализация пространства и температуры
x = np.linspace(0, L, Nx)
T = np.zeros_like(x)

# Инициализация графика
fig, ax = plt.subplots()
line, = ax.plot(x, T, color='b', label='Температура')
ax.set_xlim(0, L)
ax.set_ylim(-1.1, 0.1)

# Функция для вычисления температуры с учетом нелинейности
def T_xi(xi, t):
    """Нелинейная функция температуры"""
    return (u * beta / kappa_0) ** (1 / beta) * (xi_0 - (xi - u * t)) ** (1 / beta)

# Функция обновления графика
def update(frame):
    global T
    T_new = np.zeros_like(T)

    for i in range(1, Nx - 2):  # Обновление температуры для всех внутренних точек
        # Используем явную схему с нелинейной функцией
        xi = i / Nx  # Нормированная позиция
        T_xi_value = T_xi(xi, frame * dt)  # Вычисление значения нелинейной функции
        T_new[i] = T[i] + T_xi_value * (T[i + 1] - 2 * T[i] + T[i - 1])  # Явная схема
    
    # Граничные условия
    T_new[0] = -1  # T(0, t) = -1
    T_new[-1] = 0  # T(1, t) = 0

    T[:] = T_new[:]  # Обновление массива температуры
    
    # Обновление графика
    line.set_ydata(T)
    return line,

# Анимация
ani = FuncAnimation(fig, update, frames=Nt, blit=True, interval=50)

# Настройки графика
plt.xlabel('x')
plt.ylabel('Температура')
plt.title('Нелинейное уравнение теплопроводности')
plt.legend()

plt.show()