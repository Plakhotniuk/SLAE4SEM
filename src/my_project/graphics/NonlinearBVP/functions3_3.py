import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_nonlin_3_3.txt')


ax.plot(data[0, :], data[1, :], '.', label=f'Численное решение, 100 узлов', markersize=4)


# ax.plot(data[2, :], data[0, :], '-', label='Аналитическое решение\n y = 2th(x)')
plt.legend(fontsize=16)
ax.grid()
plt.title('y\'\' + 1/2 y y\' = -1/2y\ny(0) = 1, y(5) = 6')
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
# plt.savefig('graphic_nonlin3_yx_3.png')
plt.show()