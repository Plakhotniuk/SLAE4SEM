import matplotlib.pyplot as plt
import numpy as np

plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_nonlin_3_2.txt')

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
plt.title('Зависимость модуля максимальной ошибки от количества узлов\n'
          'y\'\' + $y\'^2$ = 0')
plt.legend()
plt.xlabel('n', fontsize=14)
plt.ylabel('max abs  error', fontsize=14)
# plt.savefig('graphic_nonlin3_max_err.png')
plt.show()