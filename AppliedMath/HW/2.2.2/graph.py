# Импортируем пакеты
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

h = 100
T = 100

x = np.linspace(0, 1, h)
t = np.linspace(0, 1, T)
y = x.copy()
xgrid, ygrid = np.meshgrid(x, y)

z0 = np.cos(np.pi * xgrid) * np.sin(5 * np.pi * ygrid) * np.exp(- 50 * pow(np.pi, 2) * 0.0001)

def z(t):
    return np.cos(np.pi * xgrid) * np.sin(5 * np.pi * ygrid) * np.exp(- 50 * pow(np.pi, 2) * 0.0001 * t)

with open('C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/f100.txt', 'r') as file:
    l = [[(num) for num in line.split(',')] for line in file]

for i in range(h * T):
    l[i].pop()

l = np.array([[float(i) for i in j] for j in l])

delta = []

def animate_func(num):
   ax.clear()
   ax.plot_surface(xgrid, ygrid, l[0 + h * num:h + h * num].T)
   ax.plot_surface(xgrid, ygrid, z(num / T))
   ax.set_title('\nTime = ' + str(np.round(t[num], decimals=2)) + ' sec')

   ax.set_xlabel('x')
   ax.set_ylabel('y')
   ax.set_zlabel('z')

# Рисуем анимацию
fig = plt.figure()
ax = plt.axes(projection='3d')
line_ani = animation.FuncAnimation(fig, animate_func, interval=1, frames=T)
plt.show()

# f = "N1000.gif"
# writergif = animation.PillowWriter(fps=100)
# line_ani.save(f, writer=writergif)

