"""Analytical continuation of Gegenbauer polynomial from Mueller/Schaefer paper."""

import numpy as np

#from .kinematics import K2, tmin, weight_BH

from scipy.special import gamma

## Original code by student Filip Bilandžija, made for 5th year
## research seminar in 2023, using algorithm from Numerical Recipes.
## Some vectorization by ChatGPT

rkh_default = 0.1  # Runge-Kutta step, increase for complicated functions

def runge_kutta4(x0, y0, x, h, aa, bb, cc, z0, dz):
    # Count number of iterations using step size or
    # step height h
    n = (int)((x - x0)/h)
    # Iterate for number of iterations
    y = y0
    for i in range(1, n + 1):
        k1 = h * hypdrv(x0, y, aa, bb, cc, z0, dz)
        k2 = h * hypdrv(x0 + 0.5 * h, y + 0.5 * k1, aa, bb, cc, z0, dz)
        k3 = h * hypdrv(x0 + 0.5 * h, y + 0.5 * k2, aa, bb, cc, z0, dz)
        k4 = h * hypdrv(x0 + h, y + k3, aa, bb, cc, z0, dz)
 
        # Update next value of y
        y += (k1 + 2 * k2 + 2 * k3 + k4)/6.
 
        # Update next value of x
        x0 += h
    return y

def hypgeo(a, b, c, z, rkh=rkh_default):
    ans, y = 0j, [0. + 0.0j, 0. + 0.0j]
    kmax=0
    if z.real * z.real + z.imag * z.imag <= 0.25:
        ans, y[1] = hypser(a, b, c, z, ans, y[1])
        return ans
    
    if z.real < 0.0:    
        z0 = -0.5 + 0.0j
    elif z.real <= 1.0:
        z0 = 0.5 + 0.0j
    else:
        z0 = -0.5j + complex(0, int(z.imag>0.0))
        
    aa, bb, cc = a, b, c # Load the global variables to pass parameters “over the head” of odeint to hypdrv.

    dz = z - z0
    y[0], y[1] = hypser(aa, bb, cc, z0, y[0], y[1])
    # odeint(yy,4,0.0,1.0,EPS,0.1,0.0001,&nok,&nbad,hypdrv,bsstep)
    yy = runge_kutta4(0.0, np.array([y[0], y[1]]), 1.0, rkh, aa, bb, cc, z0, dz)
    # The arguments to odeint are the vector of independent variables, its length, the starting
    # and ending values of the dependent variable, the accuracy parameter, an initial guess for
    # stepsize, a minimum stepsize, the (returned) number of good and bad steps taken, and the
    # names of the derivative routine and the (here Bulirsch-Stoer) stepping routine.
    return yy[0]

def hypser(a, b, c, z, series, deriv):
# Returns the hypergeometric series 2F1 and its derivative, iterating to machine accuracy. For
# |z| ≤ 1/2 convergence is quite rapid.
    deriv = 0.0 + 0.0j
    fac = complex(1.0, 0.0)
    aa, bb, cc, temp = a, b, c, fac

    for n in range(1, 1001): 
        fac = fac * (aa * bb) /cc
        deriv += fac
        fac *= 1.0/n * z
        series = temp + fac
        if series.real == temp.real and series.imag == temp.imag:
            return series, deriv
        temp = series
        aa = aa + complex(1, 0)
        bb = bb + complex(1, 0)
        cc = cc + complex(1, 0)
    
    return ArithmeticError("convergence failure in hypser")

def hypdrv(s, y, aa, bb, cc, z0, dz):
    dyds = [0.0j, 0.0j]
# Computes derivatives for the hypergeometric equation, see text equation (5.14.4).
    z = z0 + s * dz
    dyds[0] = y[1] * dz
    dyds[1] = ((aa * bb * y[0]) - (cc - (aa + bb + complex(1, 0)) * z) * y[1]) * (dz /(z * (complex(1, 0) - z)))
    return np.array([dyds[0], dyds[1]])

def heaviside_theta(cond):
    if cond > 0.:
        return True
    return False

def heaviside_theta_eq(cond):
    if cond >= 0.:
        return True
    return False


def P_j(j, x, lam, rkh=rkh_default):
    fac1 = np.power(2., j + lam - 0.5) * gamma(1 + lam + j)
    fac2 = gamma(0.5) * gamma(1. + j) * gamma(lam+0.5)
    return fac1/fac2 * (1. + x)**(lam-0.5) * hypgeo(-j - lam + 0.5, j + lam + 0.5,
                                                    lam + 0.5, (1. + x)/2., rkh)

def Q_j(j, x, lam, rkh=rkh_default):
    fac1 = np.sin(np.pi * j)/np.pi * np.power(x, -j - 1.)
    return -fac1 * hypgeo((j + 1.)/2., (j + 2.)/2., lam + 1 + j, 1 / (x ** 2), rkh)

def p_j1D(j, x, eta, lam, rkh=rkh_default):
    # Handle case where x is an array
    if isinstance(x, (list, tuple, np.ndarray)):
        res = np.zeros((len(x), len(j)), dtype=complex)

    # Precompute power matrix (eta^(-j-1)) for all j
    power_matrix = np.power(eta, -j - 1).reshape(-1, 1)  # Shape (len(j), 1)

    # Initialize result matrices
    p_j_matrix = np.zeros((len(x), len(j)), dtype=complex)
    q_j_matrix = np.zeros((len(x), len(j)), dtype=complex)

    # Compute p_j and q_j matrices
    for i, j_val in enumerate(j):
        for k, x_val in enumerate(x):
            eta_ratio = x_val / eta
            if heaviside_theta(eta - np.abs(x_val)):
                p_j_matrix[k, i] = P_j(j_val, eta_ratio, lam, rkh)
            if heaviside_theta_eq(x_val - eta):
                q_j_matrix[k, i] = Q_j(j_val, eta_ratio, lam, rkh)

    # Compute the result using broadcasting
    res += power_matrix.T * p_j_matrix
    res += power_matrix.T * q_j_matrix

    return res

def p_j2D(j, x, eta, lam, rkh=rkh_default):
    # Initialize result matrix
    res = np.zeros((len(x), len(eta), len(j)), dtype=complex)

    # Iterate over j, x, and eta to fill the result matrix
    for i, j_val in enumerate(j):
        for k, x_val in enumerate(x):
            for l, eta_val in enumerate(eta):
                eta_ratio = x_val / eta_val
                power_term = np.power(eta_val, -j_val - 1)

                if heaviside_theta(eta_val - np.abs(x_val)):
                    res[k, l, i] += power_term * P_j(j_val, eta_ratio, lam, rkh)
                
                if heaviside_theta_eq(x_val - eta_val):
                    res[k, l, i] += power_term * Q_j(j_val, eta_ratio, lam, rkh)

    return res
    
def ensure_array(x):
    if np.isscalar(x):
        return np.asarray([x])
    return np.asarray(x)


def p_j_fb(j, x, eta, lam, rkh=rkh_default):
    j = ensure_array(j)
    x = ensure_array(x)
    # Handle 2D case when eta is an array, else handle 1D case
    if isinstance(eta, (list, tuple, np.ndarray)):
        return np.moveaxis(p_j2D(j, x, eta, lam, rkh), -1, 0)
    return np.moveaxis(p_j1D(j, x, eta, lam, rkh), -1, 0)

p_j = p_j_fb
