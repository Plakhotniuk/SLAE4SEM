import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_beam1.txt')
#theta: data[1,:]
#l: data[0,:]
# ax.plot(data[0, :], data[2, :], '-', label='Аналитическое решение\n y(x) = ln(1+x) + 1')

for i in range(0, data.shape[0], 2):
    ax.plot(data[i,:], data[i+1,:], 'b.', markersize=1)

# plt.legend(fontsize=16)
ax.grid()
plt.axis('equal')
plt.title(r'$\theta$\'\' + $q^2$sin($\theta$) = 0\n$\theta$({data[0,0]}) = {data[1,0]}, $\theta$({data[0,-1]}) = {data[1,-1]}')
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)

# plt.savefig('beam1.png')
plt.show()