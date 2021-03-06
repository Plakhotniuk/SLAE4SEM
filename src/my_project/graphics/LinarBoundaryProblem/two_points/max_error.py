import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots()
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_2_err3.txt')
print(data)
x = data[:, 0]
y = data[:, 1]
x = [np.log(i) for i in x]
y = [np.log(i) for i in y]
ax.plot(x, y, '.')
t = np.polyfit(x, y, 1)
f = np.poly1d(t)
xp = np.linspace(np.min(x), np.max(x), 100)
print('---')
print(f)
plt.plot(xp, f(xp), '-', label=f'Коэффициент наклона: k = {round(f[1], 2)}')
ax.grid()
plt.legend()
plt.title('Трехдиагональная прогонка')
plt.xlabel('Логарифм количества узлов')
plt.ylabel('Логарифм максимального отклонения')
plt.savefig('graphic_max_error3.png')
plt.show()