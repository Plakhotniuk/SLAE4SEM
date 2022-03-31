import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/SLAE4SEM/cmake-build-debug/tests/test_1.txt')
# ax.plot(data[0, :], data[1, :], '-')
x = np.linspace(data[0, 0], data[0, -1], 100)
y = [np.cos(i) for i in x]
ax.plot(x, y, '-', label='Аналитическое решение' )
ax.plot(data[0, :], data[2, :], '.', label=f'Метод прогонки,\nКоличество узлов {data.shape[1]}', markersize=8)
plt.legend(fontsize=16)
ax.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
# plt.savefig('graphic_time.png')
plt.show()

