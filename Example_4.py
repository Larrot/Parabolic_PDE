import numpy as np
from Nonlinear_Parabolic1D import Nonlinear_parabolic1D


def Diff(u):
    D = 0.5*u**2
    return D


def Q(u):
    Q = 0
    return Q


def u_init(x):
    if x <= 0.5:
        u_init = (0.5-x)*40.0**(0.5)
    else:
        u_init = 0
    return u_init


def u_left(x, t):
    u_left = 1/(0.9-8*t)**(0.5)
    return u_left


def u_right(x, t):
    u_right = 0
    return u_right


L = 0.5
epsilon = 10**(-3)


N = 20

t_end = 0.111

tau = 10**(-5)


u = Nonlinear_parabolic1D(Diff, Q, u_init, u_left,
                          u_right, t_end, L, tau, N, epsilon)
np.savetxt(f'Nonlinear1_T_{t_end}_N_{N}.txt', u, fmt='%.3e', delimiter=' ')
