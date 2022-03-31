import sympy

h, A, B, C, D, E = sympy.symbols('h, A, B, C, D, E')
#second derivative
# 1 < i < n - 2
data1 = sympy.linsolve([A + B + C + D + E,
                -2*A - B + D + 2*E,
                2*A + 1/2*B + 1/2*D + 2*E - 1/h**2,
                -4/3* A - 1/6 * B + 1/6 * D + 4/3 * E,
                2/3*A + 1/24*B + 1/24*D + 2/3*E], [A, B, C, D, E])
# i = 1
data2 = sympy.linsolve([A + B + C + D + E,
                -A + C + 2*D + 3*E,
                1/2*A + 1/2*C + 2*D + 9/2*E - 1/h**2,
                -1/6* A + 1/6 * C + 4/3 * D + 9/2 * E,
                1/24*A + 1/24*C + 2/3*D + 27/8*E], [A, B, C, D, E])
#first derivative
# i = 1
data3 = sympy.linsolve([A + B + C + D + E,
                -A + C + 2*D + 3*E - 1/h,
                1/2*A + 1/2*C + 2*D + 9/2*E,
                -1/6* A + 1/6 * C + 4/3 * D + 9/2 * E,
                1/24*A + 1/24*C + 2/3*D + 27/8*E], [A, B, C, D, E])
#first derivative
# i = n - 2
data4 = sympy.linsolve([A + B + C + D + E,
                A - C - 2*D - 3*E - 1/h,
                1/2*A + 1/2*C + 2*D + 9/2*E,
                1/6* A - 1/6 * C - 4/3 * D - 9/2 * E,
                1/24*A + 1/24*C + 2/3*D + 27/8*E], [A, B, C, D, E])

#Second derivative
# i = n - 2
data5 = sympy.linsolve([A + B + C + D + E,
                A - C - 2*D - 3*E,
                1/2*A + 1/2*C + 2*D + 9/2*E - 1/h**2,
                1/6* A - 1/6 * C - 4/3 * D - 9/2 * E,
                1/24*A + 1/24*C + 2/3*D + 27/8*E], [A, B, C, D, E])

sympy.pprint(data4)
#0.91(6) = 11/12
#1.(6) = 5/3
#0.083 = 1/12
#0.83 = 5/6
