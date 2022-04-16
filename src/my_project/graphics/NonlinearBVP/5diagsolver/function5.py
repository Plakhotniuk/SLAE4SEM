import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_nonlin_5_3.txt')

# ax.plot(data[0, :], data[2, :], '-', label='Аналитическое решение\n y(x) = ln(1+x) + 1')

ax.plot(data[0, :], data[1, :], '.', label=f'Численное решение, 200 узлов', markersize=3)

plt.legend(fontsize=16)
ax.grid()
plt.axis('equal')
plt.title(f'y\'\' + (1+$y\'^2$)/(2y) = 0\ny({data[0,0]}) = {data[1,0]}, y({data[0,-1]}) = {data[1,-1]}')
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('brachistochrone5.png')
plt.show()