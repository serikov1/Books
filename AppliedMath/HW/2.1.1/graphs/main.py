import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
AutoMinorLocator)
import numpy as np
# import pandas as pd
#библиотеки

#реализация выбора
print("1 - Вывести график на экран")
print("2 - Сохранить график на рабочий стол")
proverca = int(input())

# Ua_4 = [3.33, 5.66, 6.65, 7.31, 8.40, 9.69, 10.68, 11.87, 12.92, 1377, 14.82, 15.79, 21.9, 22.30, 23.8, 25.71, 27.75, 30.13, 32.03, 34.85, 37.51, 40.18, 44.63, 47.39, 50.23, 54.06, 57.71, 60.65, 64.20, 69.28, 72.29, 77.25, 80.19]
# Ik_4 = [18.5, 34, 41, 45.8, 53.2, 61.5, 67.7, 74.9, 81.2, 86.4, 92.7, 98.3, 105.1, 86.6, 64.9, 79.7, 100.8, 125, 144.6, 171.2, 181.5, 172, 155.1, 161.5, 173.9, 196.3, 217.6, 229, 233.5, 241.5, 248.5, 265.1, 276.3]


Ik_4 = [
150, 149.74, 149.48, 149.219, 148.959, 148.699, 148.439, 148.179, 147.919, 147.659, 147.399, 147.139, 146.88, 146.62, 146.361, 146.102, 145.842, 145.583, 145.324, 145.066, 144.807, 144.549, 144.291, 144.033, 143.775, 143.518, 143.26, 143.003, 142.747, 142.49, 142.234, 141.978, 141.722, 141.467, 141.212, 140.957, 140.702, 140.448, 140.195, 139.941, 139.688, 139.435, 139.183, 138.931, 138.679, 138.428, 138.177, 137.927, 137.677, 137.428, 137.178, 136.93, 136.682, 136.434, 136.187, 135.94, 135.694, 135.448, 135.202, 134.958, 134.713, 134.47, 134.226, 133.984, 133.741, 133.5, 133.259, 133.018, 132.778, 132.539, 132.3, 132.062, 131.824, 131.587, 131.351, 131.115, 130.88, 130.645, 130.411, 130.178, 129.945, 129.713, 129.482, 129.251, 129.021, 128.792, 128.563, 128.335, 128.107, 127.881, 127.655, 127.429, 127.204, 126.981, 126.757, 126.535, 126.313, 126.092, 125.871, 125.651, 125.432, 125.214, 124.997, 124.78, 124.564, 124.348, 124.134, 123.92, 123.706, 123.494, 123.282, 123.071, 122.861, 122.652, 122.443, 122.235, 122.028, 121.821, 121.616, 121.411, 121.207, 121.003, 120.801, 120.599, 120.398, 120.197, 119.998, 119.799, 119.601, 119.403, 119.207, 119.011, 118.816, 118.621, 118.428, 118.235, 118.043, 117.851, 117.661, 117.471, 117.282, 117.093, 116.906, 116.719, 116.533, 116.347, 116.163, 115.979, 115.795, 115.613, 115.431, 115.25, 115.07, 114.89, 114.711, 114.533, 114.355, 114.179, 114.003, 113.827, 113.652, 113.478, 113.305, 113.132, 112.96, 112.789, 112.618, 112.448, 112.278, 112.11, 111.942, 111.774, 111.607, 111.441, 111.275, 111.11, 110.946, 110.782, 110.619, 110.456, 110.294, 110.133, 109.972, 109.812, 109.652, 109.493, 109.334, 109.176, 109.018, 108.861, 108.705, 108.549, 108.393, 108.238, 108.084, 107.93, 107.776, 107.623, 107.47, 107.318, 107.167, 107.015, 106.865, 106.714, 106.564, 106.415, 106.265, 106.117, 105.968, 105.82, 105.672, 105.525, 105.378, 105.232, 105.085, 104.939, 104.794, 104.648, 104.503, 104.359, 104.214, 104.07, 103.926, 103.783, 103.639, 103.496, 103.353, 103.21, 103.068, 102.926, 102.784, 102.642, 102.5, 102.358, 102.217, 102.076, 101.935, 101.794, 101.653, 101.512, 101.371, 101.231, 101.09, 100.95, 100.81, 100.669, 100.529, 100.389, 100.249, 100.108, 99.9683, 99.8281, 99.6879, 99.5477, 99.4075, 99.2672, 99.1269, 98.9866, 98.8461, 98.7056, 98.565, 98.4243, 98.2835, 98.1426, 98.0015, 97.8604, 97.719, 97.5776, 97.4359, 97.2941, 97.1521, 97.01, 96.8676, 96.725, 96.5822, 96.4392, 96.2959, 96.1525, 96.0087, 95.8647, 95.7204, 95.5758, 95.431, 95.2858, 95.1404, 94.9946, 94.8485, 94.7021, 94.5554, 94.4083, 94.2608, 94.113, 93.9648, 93.8162, 93.6672, 93.5178, 93.3681, 93.2179, 93.0673, 92.9162, 92.7647, 92.6128, 92.4604, 92.3076, 92.1543, 92.0005, 91.8462, 91.6915, 91.5362, 91.3805, 91.2242, 91.0674, 90.9101, 90.7523, 90.5939, 90.4349, 90.2755, 90.1154, 89.9548, 89.7936, 89.6319, 89.4695, 89.3066, 89.143, 88.9789, 88.8141, 88.6488, 88.4828, 88.3162, 88.149, 87.9811, 87.8126, 87.6434, 87.4736, 87.3031, 87.132, 86.9602, 86.7877, 86.6146, 86.4408, 86.2663, 86.0911, 85.9152, 85.7386, 85.5613, 85.3833, 85.2046, 85.0252, 84.8451, 84.6642, 84.4827, 84.3004, 84.1174, 83.9336, 83.7491, 83.5639, 83.3779, 83.1912, 83.0038, 82.8156, 82.6267, 82.437, 82.2465, 82.0553, 81.8634, 81.6707, 81.4772, 81.283, 81.088, 80.8923, 80.6957, 80.4985, 80.3004, 80.1016, 79.9021, 79.7017, 79.5006, 79.2988, 79.0961, 78.8927, 78.6886, 78.4837, 78.278, 78.0715, 77.8643, 77.6563, 77.4476, 77.2381, 77.0278, 76.8168, 76.605, 76.3925, 76.1792, 75.9652, 75.7504, 75.5349, 75.3186, 75.1016, 74.8839, 74.6654, 74.4461, 74.2262, 74.0055, 73.7841, 73.562, 73.3391, 73.1156, 72.8913, 72.6663, 72.4406, 72.2142, 71.9871, 71.7593, 71.5309, 71.3017, 71.0719, 70.8414, 70.6102, 70.3783, 70.1458, 69.9127, 69.6788, 69.4444, 69.2093, 68.9735, 68.7372, 68.5002, 68.2626, 68.0243, 67.7855, 67.5461, 67.306, 67.0654, 66.8242, 66.5824, 66.3401, 66.0972, 65.8537, 65.6097, 65.3651, 65.12, 64.8744, 64.6282, 64.3816, 64.1344, 63.8867, 63.6386, 63.3899, 63.1408, 62.8912, 62.6411, 62.3906, 62.1396, 61.8882, 61.6364, 61.3841, 61.1315, 60.8784, 60.6249, 60.3711, 60.1168, 59.8622, 59.6072, 59.3519, 59.0962, 58.8402, 58.5838, 58.3272, 58.0702, 57.8129, 57.5553, 57.2975, 57.0393, 56.7809, 56.5223, 56.2634, 56.0042, 55.7448, 55.4852, 55.2254, 54.9654, 54.7052, 54.4449, 54.1843, 53.9236, 53.6627, 53.4017, 53.1406, 52.8793, 52.6179, 52.3565, 52.0949, 51.8332, 51.5715, 51.3097, 51.0478, 50.7859, 50.524, 50.262, 50,
]
Ua_4 = [x for x in range(len(Ik_4))]
# x1 = 2.51, 2.58, 2.67, 2.78, 2.87, 2.94, 3.03, 3.12, 3.23, 3.32, 3.41, 3.48, 3.57, 3.68, 3.75, 3.86, 3.93, 4.02, 4.13, 4.22, 4.31, 4.38, 4.47, 4.56, 4.67, 4.74, 4.85, 4.92, 5.01, 5.12, 5.21
# y1 = 50.2, 50.97, 51.75, 55.15, 56.17, 52.81, 56.08, 56.95, 58.06, 59.39, 60.32, 63.7, 60.89, 66.29, 65.13, 66.79, 65.88, 71.38, 68.32, 71.89, 75.27, 76.48, 76.32, 75.79, 76.88, 82.41, 81.96, 83.55, 87.16, 89.21, 88.89

#x2 = 206.06928,	397.01424,	579.65056,	745.15392,	884.7,	989.46448,	1044.576039
#y2 = 4.5, 9, 16.5, 21, 27, 28.5, 30

#x3 = 206.06928,	397.01424,	579.65056,	745.15392,	884.7,	989.46448,	1044.576039
#y3 = 9, 16.5, 27, 34.5, 40.5, 45, 46.5

#x4 = 206.06928,	397.01424,	579.65056,	745.15392,	884.7,	989.46448,	1044.576039
#y4 = 10.5, 21, 33, 45, 52.5, 57, 60

#x5 = 206.06928,	397.01424,	579.65056,	745.15392,	884.7,	989.46448,	1044.576039
#y5 = 12, 28.5, 43.5, 57, 66, 73.5, 78

#x6 = 206.06928,	397.01424,	579.65056,	745.15392,	884.7,	989.46448,	1044.576039
#y6 = 13.5, 33, 49.5, 66, 78, 88.5, 91.5



#xerr1=np.array([3.09104, 5.95521, 8.69476, 11.1773, 13.2705, 14.842])
#yerr1=np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5])

#xerr2=np.array([3.09104, 5.95521, 8.69476, 11.1773, 13.2705, 14.842, 15.6686])
#yerr2=np.array([1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5])

#xerr3=np.array([3.09104, 5.95521, 8.69476, 11.1773, 13.2705, 14.842, 15.6686])
#yerr3=np.array([1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5])

#xerr4=np.array([3.09104, 5.95521, 8.69476, 11.1773, 13.2705, 14.842, 15.6686])
#yerr4=np.array([1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5])

#xerr5=np.array([3.09104, 5.95521, 8.69476, 11.1773, 13.2705, 14.842, 15.6686])
#yerr5=np.array([1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5])

#xerr6=np.array([3.09104, 5.95521, 8.69476, 11.1773, 13.2705, 14.842, 15.6686])
#yerr6=np.array([1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5])




tochki1 = np.linspace(Ua_4[0], Ua_4[-1], 10000)
#tochki2 = np.linspace(x2[0], x2[-1], 10000)
#tochki3 = np.linspace(x3[0], x3[-1], 10000)
#tochki4 = np.linspace(x4[0], x4[-1], 10000)
#tochki5 = np.linspace(x5[0], x5[-1], 10000)
#tochki6 = np.linspace(x6[0], x6[-1], 10000)


z1 = np.polyfit(Ua_4, Ik_4, 2)
p1 = np.poly1d(z1)

#z2 = np.polyfit(x2, y2, 1)
#p2 = np.poly1d(z2)

#z3 = np.polyfit(x3, y3, 1)
#p3 = np.poly1d(z3)

#z4 = np.polyfit(x4, y4, 1)
#p4 = np.poly1d(z4)

#z5 = np.polyfit(x5, y5, 1)
#p5 = np.poly1d(z5)

#z6 = np.polyfit(x6, y6, 1)
#p6 = np.poly1d(z6)




#fig, ax = plt.subplots(figsize=(10, 7))
if proverca == 1:
    fig, ax = plt.subplots(figsize=(10, 7))
else:
    fig, ax = plt.subplots(figsize=(10, 7), dpi = 600)
#для корректоного вывода на экан, не трогать


#plt.axis([16,32,0,4.6])
##обрезка координат


ax.set_title("P for t = 10 days", fontsize=16)                    #название графика
#названия и имена



#ax.set_yscale('log')
#логарифмический масштаб для оси Y



ax.grid(which="major", linewidth=1.3)                               #мажорная сетка
ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5) #минорная сетка
#создаём сетку для графика


ax.plot(Ua_4, Ik_4,"r.", markersize=8, label = 'Ток на образце 0.2 А')
#ax.plot(x2, y2,"b.", markersize=8, label = 'Ток на образце 0.4 А' )
#ax.plot(x3, y3,"g.", markersize=8, label = 'Ток на образце 0.6 А' )
#ax.plot(x4, y4,"y.", markersize=8, label = 'Ток на образце 0.8 А')
#ax.plot(x5, y5,"k.", markersize=8, label = 'Ток на образце 1.0 А' )
#ax.plot(x6, y6,"m.", markersize=8, label = 'Ток на образце 1.2 А' )
#строительство графика на рисунке
# ax.plot(tochki1, p1(tochki1), 'b--', label = '')
#ax.plot(tochki2, p2(tochki2), 'b--', label = '')
#ax.plot(tochki3, p3(tochki3), 'g--', label = '')
#ax.plot(tochki4, p4(tochki4), 'y--', label = '')
#ax.plot(tochki5, p5(tochki5), 'k--', label = '')
#ax.plot(tochki6, p6(tochki6), 'm--', label = '')


#ax.plot(x2, p2(x2), 'g--', label = 'Максимальный наклон кривой')
#ax.plot(x3, p3(x3), 'b--')
#в скобказ указываем точки графика для которого сторим линию тренда и функцию полинома
#строительство линии тренда
#ax.set_xticks(numpy.arange(0, 1000, 10))
#ax.set_yticks(numpy.arange(0, 1., 0.1))
#ax.legend()
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='major', length=10, width=1)
ax.tick_params(which='minor', length=5, width=1)
#строительство минорной и мажорной сетки

#ax.errorbar(x1, y1,
#xerr=xerr1,
#yerr=yerr1,
#fmt='.', color='red', markersize=5)

#ax.errorbar(x2, y2,
#xerr=xerr2,
#yerr=yerr2,
#fmt='.', color='blue', markersize=5)

#ax.errorbar(x3, y3,
#xerr=xerr3,
#yerr=yerr3,
#fmt='.', color='green', markersize=5)

#ax.errorbar(x4, y4,
#xerr=xerr4,
#yerr=yerr4,
#fmt='.', color='yellow', markersize=5)

#ax.errorbar(x5, y5,
#xerr=xerr5,
#yerr=yerr5,
#fmt='.', color='black', markersize=5)

#ax.errorbar(x6, y6,
#xerr=xerr6,
#yerr=yerr6,
#fmt='.', color='purple', markersize=5)






if proverca == 1:
    plt.show()
else:
    plt.savefig("10days.png")
