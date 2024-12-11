import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Параметры задачи
Lx, Ly, Lz = 10.0, 10.0, 10.0  # размеры бруска (м)
Nx, Ny, Nz = 20, 20, 20         # количество ячеек по осям x, y, z
Q = 1.0                         # плотность теплового потока (Вт/м^3)
k = 1.0                         # теплопроводность (Вт/(м·К))

# Размеры ячеек
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)
dz = Lz / (Nz - 1)

# Инициализация температуры (начальное приближение)
T = np.zeros((Nx, Ny, Nz))

# Граничные условия (например, на границах температура 0)
T[:, 0, :] = 0.0  # y=0
T[:, -1, :] = 0.0  # y=L
T[0, :, :] = 0.0   # x=0
T[-1, :, :] = 0.0  # x=L
T[:, :, 0] = 0.0   # z=0
T[:, :, -1] = 0.0  # z=L

# Решение уравнения Пуассона методом конечных разностей
def solve_poisson(T, Q, k, dx, dy, dz, tol=1e-6, max_iter=10000):
    frames = []  # Список для хранения кадров для анимации
    for _ in range(max_iter):
        T_new = np.copy(T)
        
        # Обновление температуры для всех внутренних точек
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                for l in range(1, Nz-1):
                    T_new[i, j, l] = (T[i+1, j, l] + T[i-1, j, l]) * dy**2 * dz**2 + \
                                     (T[i, j+1, l] + T[i, j-1, l]) * dx**2 * dz**2 + \
                                     (T[i, j, l+1] + T[i, j, l-1]) * dx**2 * dy**2
                    T_new[i, j, l] /= 2 * (dx**2 * dy**2 + dy**2 * dz**2 + dz**2 * dx**2)
                    T_new[i, j, l] += Q / k * (dx**2 * dy**2 * dz**2)
        
        # Добавляем текущий кадр (сечение по оси X, Y и Z)
        frames.append([T_new[:, Ny//2, :], T_new[:, :, Nz//2], T_new[Nx//2, :, :]])
        
        # Проверка на сходимость
        if np.max(np.abs(T_new - T)) < tol:
            break
        
        T = np.copy(T_new)
    
    return frames

# Решение задачи и получение кадров для анимации
frames = solve_poisson(T, Q, k, dx, dy, dz)

# Настройка графика для анимации
fig, ax = plt.subplots(1, 3, figsize=(15, 5))

# Функция обновления кадров для анимации
def update(frame):
    # Обновление каждого сечения (по осям X, Y, Z)
    im1 = ax[0].imshow(frame[0], cmap='plasma', origin='lower', extent=[0, Lx, 0, Lz])
    ax[0].set_title("Сечение по оси X")
    ax[0].set_xlabel("x (м)")
    ax[0].set_ylabel("z (м)")
    
    im2 = ax[1].imshow(frame[1], cmap='plasma', origin='lower', extent=[0, Ly, 0, Lz])
    ax[1].set_title("Сечение по оси Y")
    ax[1].set_xlabel("x (м)")
    ax[1].set_ylabel("z (м)")
    
    im3 = ax[2].imshow(frame[2], cmap='plasma', origin='lower', extent=[0, Lx, 0, Ly])
    ax[2].set_title("Сечение по оси Z")
    ax[2].set_xlabel("x (м)")
    ax[2].set_ylabel("y (м)")
    
    # Установка цветовых полос
    fig.colorbar(im1, ax=ax[0], label="Температура (K)")
    fig.colorbar(im2, ax=ax[1], label="Температура (K)")
    fig.colorbar(im3, ax=ax[2], label="Температура (K)")
    
    return [im1, im2, im3]

# Создание анимации
ani = animation.FuncAnimation(fig, update, frames=frames, interval=100, blit=False)

# Отображение анимации
plt.tight_layout()
plt.show()