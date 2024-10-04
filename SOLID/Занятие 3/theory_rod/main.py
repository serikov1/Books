import matplotlib.pyplot as plt
import numpy as np
from math import sinh
from math import pi
from math import sqrt

# Создаём экземпляр класса figure и добавляем к Figure область Axes
fig, ax = plt.subplots()
# Добавим заголовок графика
ax.set_title('Распределение температуры')
# Название оси X:
ax.set_xlabel('x, м')
# Название оси Y:
ax.set_ylabel('T, К')
# Начало и конец изменения значения X, разбитое на 100 точек
X = np.linspace(0, 1, 100)

# Построение
T_0 = 300
T_1 = 373
T_2 = 273
l = 1
a = 50
d = 0.01
hi = 200
c = 900
ro = 2700

p = 2 * pi * d/2
S = pi * d**2 / 4
b = sqrt((a * p) / S)
beta = b / sqrt(hi)
print(beta)
T = [T_0 + ((T_1 - T_0) * sinh(beta*(l - x)) + (T_2 - T_0) * sinh(beta * x))/(sinh(beta * l)) for x in X]
# Вывод графика
ax.plot(X, T)
plt.grid(True)
plt.savefig('T(x)-theoretical')
plt.show()
