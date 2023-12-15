import cmath
from matplotlib import pyplot
import matplotlib
import numpy as np

matplotlib.use('TkAgg')
y = [
6.14903,
4.7739,
1.93254,
8.07588,
0.607833,
6.99694,
3.74315,
2.81437,
7.62679,
0.479642,
]
x = [x * 10 for x in range(len(y))]

pyplot.scatter(x, y, s=10., color='red')
pyplot.legend(['$\\Delta$ x = 1'])
pyplot.grid()

pyplot.show()
