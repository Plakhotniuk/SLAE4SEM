import matplotlib.pyplot as plt
import numpy as np

plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_nonlin_3.txt')

n_iterations = data.shape[0] - 2

y_analytic = data[1, :]

max_abs_errors = [np.max(np.abs(data[i,:] - y_analytic)) for i in range(2, data.shape[0])]
print(max_abs_errors)
ax.plot([i for i in range(n_iterations)], max_abs_errors, '-')

ax.grid()
plt.title('Зависимость модуля максимальной ошибки от количества итераций\n'
          'y\'\' + y y\' = 0; y(0) = 0; y(1) = 2 th2')
plt.xlabel('n', fontsize=14)
plt.ylabel('max abs  error', fontsize=14)
# plt.savefig('graphic_nonlin3_max_err.png')
plt.show()