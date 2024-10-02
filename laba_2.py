import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

# Заданные параметры
L = 1.0  # длина стержня
T_max = 4.0  # конечное время
h = 0.02  # шаг по пространству
tau = 0.0001  # шаг по времени
x = np.arange(0, L + h, h)  # сетка по пространству
t = np.arange(0, T_max + tau, tau)  # сетка по времени
N_x = len(x)
N_t = len(t)

# Коэффициент alpha
alpha = tau / h**2

# Инициализация решения
T_implicit = np.zeros((N_t, N_x))

# Начальные условия: T(x, 0) = 0
T_implicit[0, :] = 0

# Функция для граничного условия на левой границе
def left_boundary(ti):
    return ti / (1 + ti)

# Коэффициенты для трехдиагональной матрицы
a = -alpha
b = 1 + 2 * alpha
c = -alpha

# Создание матрицы для метода прогонки
A = np.zeros((3, N_x))
A[0, 1:] = c  # верхняя диагональ
A[1, :] = b  # главная диагональ
A[2, :-1] = a  # нижняя диагональ

# Неявная схема
for n in range(0, N_t - 1):
    # Правая часть уравнения: предыдущий временной слой
    d = T_implicit[n, :].copy()

    # Граничные условия
    d[0] = left_boundary(t[n + 1])
    d[-1] = 0

    # Решение системы уравнений методом прогонки
    T_implicit[n + 1, :] = solve_banded((1, 1), A, d)

# Построение графиков для нескольких моментов времени
plt.figure(figsize=(10, 6))
for time_step in [0, int(1/tau), int(2/tau), int(3/tau), int(4/tau)]:
    plt.plot(x, T_implicit[time_step, :], label=f't = {time_step * tau:.1f} s')

plt.title('Решение уравнения теплопроводности - неявная схема')
plt.xlabel('x')
plt.ylabel('T(x,t)')
plt.legend()
plt.grid(True)
plt.show()