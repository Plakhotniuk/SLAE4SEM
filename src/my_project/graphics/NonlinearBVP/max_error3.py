import matplotlib.pyplot as plt
import numpy as np

plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_nonlin_3.txt')

n_iterations = data.shape[0] - 2

y_analytic = data[1, :]

max_abs_errors = [np.max(np.abs(data[i,:] - y_analytic)) for i in range(2, data.shape[0])]
print(max_abs_errors)
# ax.plot([i for i in range(n_iterations)], max_abs_errors, '-')
x = [i for i in range(1, n_iterations + 1)]

x = [np.log(i) for i in x]
y = [np.log(i) for i in max_abs_errors]
ax.plot(x, y, '.')
t = np.polyfit(x, y, 1)
f = np.poly1d(t)
xp = np.linspace(np.min(x), np.max(x), 100)
print('---')
print(f)
# plt.plot(xp, f(xp), '-', label=f'Коэффициент наклона: k = {round(f[1], 2)}')

ax.grid()
plt.title('Зависимость модуля максимальной ошибки от количества итераций\n'
          'y\'\' + y y\' = 0; y(0) = 0; y(1) = 2 th2')
plt.xlabel('n', fontsize=14)
plt.ylabel('max abs  error', fontsize=14)
# plt.savefig('graphic_nonlin3_max_err.png')
plt.show()