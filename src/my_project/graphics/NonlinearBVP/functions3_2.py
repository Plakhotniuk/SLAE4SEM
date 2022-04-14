import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_nonlin_3_2.txt')

n_iterations = data.shape[0] - 2
for i in range(2, data.shape[0]):
    ax.plot(data[0, :], data[i, :], 'x', label=f'{i - 1}', markersize=3)


ax.plot(data[0, :], data[1, :], '-', label='Аналитическое решение\n y = 2th(x)')
plt.legend(fontsize=16)
ax.grid()
plt.title('y\'\' + y y\' = 0\ny(0) = 0, y(1) = 2 th2')
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
# plt.savefig('graphic_nonlin3_yx.png')
plt.show()