import numpy as np


def Nonlinear_parabolic1D(D, Q, u_init, u_left, u_right, t_end, L, tau, N, epsilon):
    """
    Implicit scheme for numerical solution of the Dirichle 
    problem for one-dimensional parabolic equation.


    Parameters
    ----------
    D : function
        diffusion coefficient D=D(u).
    Q : function
        source term Q=Q(u).
    u_init : function
        initial condition u=u(0,x).
    u_left : function
        Left boundary condition u=u(t,0).
    u_right : function
        Right boundary condition u=u(t,L).
    t_end : float 
        End time of computations.
    L : float 
        The length of spatial area.
    tau : float 
        the length of time step.
    N : int 
        the number of nodes.
    epsilon : float 
        Accuracy of iterative process.    

    Returns function u(x, t_end)
    -------
    """
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
             Right side of equations.

         Returns solution array x
         -------
        """
        size = len(d)
        x = np.zeros(size)
        # Forward steps
        for i in range(size-1):
            w = a[i]/b[i]
            b[i+1] -= w*c[i]
            d[i+1] -= w*d[i]
        # Backward steps
        x[size-1] = d[size-1]/b[size-1]
        for i in range(size-2, -1, -1):
            x[i] = (d[i]-c[i]*x[i+1])/b[i]

        return x
    # Spatial step
    h = L/(N-1)

    # const

    C = tau/h**2
    # Spatial grid
    x = np.linspace(0, L, N)
    # Initial time
    t = 0

    # Solution function
    u = np.zeros(N)
    # Array for finding error
    error = np.zeros(N)
    # Initial condition
    for i in range(N):
        u[i] = u_init(x[i])
    # Array for iterative process
    u_iter = u.copy()
    # Boundary conditions
    # Left
    u[0] = u_left(0, t)
    # Right
    u[N-1] = u_right(0, t)

    # Array for RSE
    RP = np.zeros(N-2)
    # Arrays for tridiagonal matrix
    a = np.zeros(N-3)
    b = np.zeros(N-2)
    c = np.zeros(N-3)

    while t < t_end:
        # Iterative process
        while True:
            u_iter_old = u_iter.copy()
            # Filling RSE
            for i in range(1, N-3):
                RP[i] = u[i+1]+Q(u_iter[i+1])*tau
            RP[0] = u[1]+tau*Q(u_iter[1])+C/2 * \
                (D(u_iter[1])+D(u_iter[0]))*u_left(0, t+tau)
            RP[N-3] = u[N-2]+tau*Q(u_iter[N-2])+C/2*(D(u_iter[N-1]) +
                                                     D(u_iter[N-2]))*u_right(L, t+tau)

            # Filling lower diagonal
            for i in range(N-3):
                a[i] = -C/2*(D(u_iter[i+2])+D(u_iter[i+1]))
            # Filling middle diagonal
            for i in range(N-2):
                b[i] = 1 + C/2*(D(u_iter[i+2])+2*D(u_iter[i+1])+D(u_iter[i]))
            # Filling upper diagonal
            for i in range(N-3):
                c[i] = -C/2*(D(u_iter[i+2])+D(u_iter[i+1]))
            # Solving system of equations
            u_iter[1:N-1] = Tridiagonal_solver(a, b, c, RP)
            # Finding error
            for i in range(N):
                error[i] = abs(u_iter_old[i]-u_iter[i])
            if max(error) < epsilon:
                break
        error = np.zeros(N)
        u = u_iter.copy()
        t += tau
    return u
