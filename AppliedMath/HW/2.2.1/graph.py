import matplotlib.pyplot as plt
from collections import deque
import random
from matplotlib.animation import FuncAnimation

import numpy as np

t = np.linspace(0, 239, 239)

# Create an empty plot
fig = plt.figure()
ax = plt.axes(xlim=(0, 500), ylim=(40, 160))
line, = ax.plot([], [], lw=3)

plt.ylabel("P, атм")
plt.xlabel("x, м")

ax.grid(which="major", linewidth=1.3)  # мажорная сетка
ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5)  # минорная сетка
with open('C:/PHYSTECH/Books/AppliedMath/HW/2.2.1/preasure.txt', 'r') as file:
    l = [[(num) for num in line.split(',')] for line in file]

# print(len(l))
for i in range(239):
    l[i].pop()

l = np.array([[float(i) for i in j] for j in l])

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = np.linspace(0, 500, 500)
    y = list(l[i])
    # print(y)
    line.set_data(x, y)
    ax.set_title('\nTime = ' + str(np.round(t[i], decimals=2)) + ' hour')
    return line,


anim = FuncAnimation(fig, animate, init_func=init, frames=239,interval=100, blit=True)
anim.save('p.gif')