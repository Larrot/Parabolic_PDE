import numpy as np


def Tridiagonal_solver(a, b, c, d):
    """
     Algorithm to solve tridiagonal systems of equations

     Parameters
     ----------
     a : array
         Lower diagonal.
     b : array
         Middle diagonal.
     c : array
         Upper diagonal.
     d : array
         right side of equation.

     Returns solution array x
     -------
    """
    size = len(d)
    x = np.zeros(size)
    for i in range(size-1):
        w = a[i]/b[i]
        b[i+1] -= w*c[i]
        d[i+1] -= w*d[i]
    x[size-1] = d[size-1]/b[size-1]
    for i in range(size-2, -1, -1):
        x[i] = (d[i]-c[i]*x[i+1])/b[i]

    return x


"""
Testing Tridiagonal_solver via built-in function np.linalg.solve
"""
a = np.array([1, 2, 3.])
b = np.array([2, 3, 4., 5.])
c = np.array([2., 3., 4.])
d = np.array([3., 4., 5., 6.])

M = np.diag(a, -1)+np.diag(b)+np.diag(c, 1)

print(np.linalg.solve(M, d))
print(Tridiagonal_solver(a, b, c, d))
