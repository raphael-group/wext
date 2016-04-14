#!/usr/bin/env python

import numpy as np, numpy.linalg
from scipy.optimize import fsolve
from scipy.stats import norm

def saddlepoint_intersection_2(z, m_x, m_y, p_x, p_y, verbose=False):

    # Precompute commonly used constants.

    q_x = 1.0-p_x
    q_y = 1.0-p_y
    z += 1

    # Define K.

    def L(r, s, t):
        return np.exp(r+s+t)*p_x*p_y + np.exp(r)*p_x*q_y + np.exp(s)*q_x*p_y + q_x*q_y

    def K(r, s, t):
        return np.sum(np.log(L(r, s, t)))

    # Define K'.

    def L_r(r, s, t):
        return np.exp(r+s+t)*p_x*p_y + np.exp(r)*p_x*q_y

    def L_s(r, s, t):
        return np.exp(r+s+t)*p_x*p_y + np.exp(s)*q_x*p_y

    def L_t(r, s, t):
        return np.exp(r+s+t)*p_x*p_y

    def K_r(r, s, t):
        return np.sum(L_r(r, s, t)/L(r, s, t))

    def K_s(r, s, t):
        return np.sum(L_s(r, s, t)/L(r, s, t))

    def K_t(r, s, t):
        return np.sum(L_t(r, s, t)/L(r, s, t))

    def dK(r, s, t):
        return np.array([K_r(r, s, t), K_s(r, s, t), K_t(r, s, t)])

    # Define K''.

    def L_rr(r, s, t):
        return np.exp(r+s+t)*p_x*p_y + np.exp(r)*p_x*q_y

    def L_rs(r, s, t):
        return np.exp(r+s+t)*p_x*p_y

    def L_rt(r, s, t):
        return np.exp(r+s+t)*p_x*p_y

    def L_ss(r, s, t):
        return np.exp(r+s+t)*p_x*p_y + np.exp(s)*q_x*p_y

    def L_st(r, s, t):
        return np.exp(r+s+t)*p_x*p_y

    def L_tt(r, s, t):
        return np.exp(r+s+t)*p_x*p_y

    def K_rr(r, s, t):
        return np.sum(L_rr(r, s, t)/L(r, s, t) - (L_r(r, s, t)**2)/(L(r, s, t)**2))

    def K_rs(r, s, t):
        return np.sum(L_rs(r, s, t)/L(r, s, t) - (L_r(r, s, t)*L_s(r, s, t))/(L(r, s, t)**2))

    def K_rt(r, s, t):
        return np.sum(L_rt(r, s, t)/L(r, s, t) - (L_r(r, s, t)*L_t(r, s, t))/(L(r, s, t)**2))

    def K_ss(r, s, t):
        return np.sum(L_ss(r, s, t)/L(r, s, t) - (L_s(r, s, t)**2)/(L(r, s, t)**2))

    def K_st(r, s, t):
        return np.sum(L_st(r, s, t)/L(r, s, t) - (L_s(r, s, t)*L_t(r, s, t))/(L(r, s, t)**2))

    def K_tt(r, s, t):
        return np.sum(L_tt(r, s, t)/L(r, s, t) - (L_t(r, s, t)**2)/(L(r, s, t)**2))

    def d2K(r, s, t):
        return np.array([[K_rr(r, s, t), K_rs(r, s, t), K_rt(r, s, t)],
                         [K_rs(r, s, t), K_ss(r, s, t), K_st(r, s, t)],
                         [K_rt(r, s, t), K_st(r, s, t), K_tt(r, s, t)]])

    # Define K_x and its derivatives.

    def L_x(r):
        return np.exp(r)*p_x + q_x

    def L_xr(r):
        return np.exp(r)*p_x

    def L_xrr(r):
        return np.exp(r)*p_x

    def K_x(r):
        return np.sum(np.log(L_x(r)))

    def dK_x(r):
        return np.sum(L_xr(r)/L_x(r))

    def d2K_x(r):
        return np.sum(L_xrr(r)/L_x(r) - (L_xr(r)**2)/(L_x(r)**2))

    # Define K_y and its derivatives.

    def L_y(s):
        return np.exp(s)*p_y + q_y

    def L_ys(s):
        return np.exp(s)*p_y

    def L_yss(s):
        return np.exp(s)*p_y

    def K_y(s):
        return np.sum(np.log(L_y(s)))

    def dK_y(s):
        return np.sum(L_ys(s)/L_y(s))

    def d2K_y(s):
        return np.sum(L_yss(s)/L_y(s) - (L_ys(s)**2)/(L_y(s)**2))

    # Solve K_xr(r_h) = m_x.

    def dK_x_wrapper((r, )):
        return np.array([dK_x(r)-m_x])

    def d2K_x_wrapper((r, )):
        return np.array([d2K_x(r)])

    r_h, = fsolve(dK_x_wrapper, 1.0, fprime=d2K_x_wrapper)

    # Solve K_ys(s_h) = m_y.

    def dK_y_wrapper((s, )):
        return np.array([dK_y(s)-m_y])

    def d2K_y_wrapper((s, )):
        return np.array([d2K_y(s)])

    s_h, = fsolve(dK_y_wrapper, 1.0, fprime=d2K_y_wrapper)

    # Solve dK(r_t, s_t, t_t) = (m_x, m_y, z-0.5).

    def dK_wrapper((r, s, t)):
        return dK(r, s, t)-np.array([m_x, m_y, z-0.5])

    def d2K_wrapper((r, s, t)):
        return d2K(r, s, t)

    r_t, s_t, t_t = fsolve(dK_wrapper, np.ones(3), fprime=d2K_wrapper)

    # Compute saddlepoint approximation.

    w_t = np.sign(t_t)*np.sqrt(2.0*((K_x(r_h)+K_y(s_h))-(m_x*r_h+m_y*s_h) - (K(r_t, s_t, t_t)-(m_x*r_t+m_y*s_t+(z-0.5)*t_t))))
    u_t = 2.0*np.sinh(0.5*t_t)*np.sqrt(np.linalg.det(d2K(r_t, s_t, t_t))/(d2K_x(r_h)*d2K_y(s_h)))

    return 1-(norm.cdf(-w_t)-norm.pdf(w_t)*(1.0/w_t-1.0/u_t))
