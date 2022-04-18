import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

x = []
y = []

with open(f'/Users/arseniy/Desktop/SLAE4SEM/cmake-build-debug/tests/test_beam9.txt', 'r') as file:
    for ind, lines in enumerate(file):
        if not ind % 2:
            x.append([float(i) for i in lines.split()])
        else:
            y.append([float(i) for i in lines.split()])

fig = plt.figure()
x_lim = [min([min(i) for i in x]), max([max(i) for i in x])]
y_lim = [min([min(i) for i in y]), max([max(i) for i in y])]

ax = plt.axes()
ax.set_xlim(-1, 1)
ax.set_ylim(-0.1, 1)

line, = ax.plot([], [], color='b')
plt.grid()


def init():
    line.set_data([], [])
    return line,


def animate(i):
    plt.title(r'$\theta$'+ f'\'\'' +
              r' + $q^2$sin($\theta$) = 0, ' +
              r'$\theta$(0) = 0, ' + r'$\theta$(1) = $\pi$'
              f'\nq = (P$L^2$)/(EI) = {round(0.025*i, 2)}')
    line.set_data(x[i], y[i])
    return line,


anim = FuncAnimation(fig, animate, init_func=init,
                     frames=170, interval=150, blit=True)

anim.save('gif9.gif', writer='pillow')
