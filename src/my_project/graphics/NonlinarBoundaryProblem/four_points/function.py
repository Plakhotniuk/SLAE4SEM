import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_4_func2.txt')
# ax.plot(data[0, :], data[1, :], '-')
x = data[:, 1]
y = data[:, 2]
ax.plot(x, y, '-', label='Аналитическое решение')
ax.plot(x, data[:, 0], 'x', label=f'Метод прогонки,\nКоличество узлов {data.shape[0]}', markersize=2)
plt.legend(fontsize=16)
ax.grid()
plt.title('y = 2*x - M_PI + M_PI*cos(x)')
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('graphic_yx2.png')
plt.show()
