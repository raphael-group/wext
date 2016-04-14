#!/usr/bin/env python

import numpy as np, numpy.linalg
from scipy.optimize import fsolve
from scipy.stats import norm

def saddlepoint_exclusive_2(z, m_x, m_y, p_x, p_y, verbose=False):

    # Precompute commonly used constants.

    q_x = 1.0-p_x
    q_y = 1.0-p_y
    m = np.array([m_x, m_y, z-0.5])

    # Define K.

    def L(r, s, t):
        return np.exp(r+t)*p_x*q_y + np.exp(s+t)*q_x*p_y + np.exp(r+s)*p_x*p_y + q_x*q_y

    def K(r, s, t):
        return np.sum(np.log(L(r, s, t)))

    # Define K'.

#    def L_r(r, s, t):
#        return np.exp(r+t)*p_x*q_y + np.exp(r+s)*p_x*p_y

#    def L_s(r, s, t):
#        return np.exp(s+t)*q_x*p_y + np.exp(r+s)*p_x*p_y

#    def L_t(r, s, t):
#        return np.exp(r+t)*p_x*q_y + np.exp(s+t)*q_x*p_y

#    def K_r(r, s, t):
#        return np.sum(L_r(r, s, t)/L(r, s, t))

#    def K_s(r, s, t):
#        return np.sum(L_s(r, s, t)/L(r, s, t))

#    def K_t(r, s, t):
#        return np.sum(L_t(r, s, t)/L(r, s, t))

#    def dK(r, s, t):
#        return np.array([K_r(r, s, t), K_s(r, s, t), K_t(r, s, t)])

    def dK_(r, s, t):

        a = np.exp(r+t)*p_x*q_y
        b = np.exp(s+t)*q_x*p_y
        c = np.exp(r+s)*p_x*p_y
        d = q_x*q_y

        L_i_ = 1.0/(a+b+c+d)
        L_r_ = a + c
        L_s_ = b + c
        L_t_ = a + b

        K_r_ = np.sum(L_r_*L_i_)
        K_s_ = np.sum(L_s_*L_i_)
        K_t_ = np.sum(L_t_*L_i_)

        return np.array([K_r_, K_s_, K_t_])

    # Define K''.

#    def L_rr(r, s, t):
#        return np.exp(r+t)*p_x*q_y + np.exp(r+s)*p_x*p_y

#    def L_rs(r, s, t):
#        return np.exp(r+s)*p_x*p_y

#    def L_rt(r, s, t):
#        return np.exp(r+t)*p_x*q_y

#    def L_ss(r, s, t):
#        return np.exp(s+t)*q_x*p_y + np.exp(r+s)*p_x*p_y

#    def L_st(r, s, t):
#        return np.exp(s+t)*q_x*p_y

#    def L_tt(r, s, t):
#        return np.exp(r+t)*p_x*q_y + np.exp(s+t)*q_x*p_y

#    def K_rr(r, s, t):
#        return np.sum(L_rr(r, s, t)/L(r, s, t) - (L_r(r, s, t)**2)/(L(r, s, t)**2))

#    def K_rs(r, s, t):
#        return np.sum(L_rs(r, s, t)/L(r, s, t) - (L_r(r, s, t)*L_s(r, s, t))/(L(r, s, t)**2))

#    def K_rt(r, s, t):
#        return np.sum(L_rt(r, s, t)/L(r, s, t) - (L_r(r, s, t)*L_t(r, s, t))/(L(r, s, t)**2))

#    def K_ss(r, s, t):
#        return np.sum(L_ss(r, s, t)/L(r, s, t) - (L_s(r, s, t)**2)/(L(r, s, t)**2))

#    def K_st(r, s, t):
#        return np.sum(L_st(r, s, t)/L(r, s, t) - (L_s(r, s, t)*L_t(r, s, t))/(L(r, s, t)**2))

#    def K_tt(r, s, t):
#        return np.sum(L_tt(r, s, t)/L(r, s, t) - (L_t(r, s, t)**2)/(L(r, s, t)**2))

#    def d2K(r, s, t):
#        return np.array([[K_rr(r, s, t), K_rs(r, s, t), K_rt(r, s, t)],
#                         [K_rs(r, s, t), K_ss(r, s, t), K_st(r, s, t)],
#                         [K_rt(r, s, t), K_st(r, s, t), K_tt(r, s, t)]])

    def d2K_(r, s, t):

        a = np.exp(r+t)*p_x*q_y
        b = np.exp(s+t)*q_x*p_y
        c = np.exp(r+s)*p_x*p_y
        d = q_x*q_y

        L_i_ = 1.0/(a+b+c+d)
        L_r_ = a + c
        L_s_ = b + c
        L_t_ = a + b

        L_ii_ = L_i_**2
        L_rr_ = L_r_
        L_rs_ = c
        L_rt_ = a
        L_ss_ = L_s_
        L_st_ = b
        L_tt_ = L_t_

        K_rr_ = np.sum(L_rr_*L_i_ - L_r_**2*L_ii_)
        K_rs_ = np.sum(L_rs_*L_i_ - L_r_*L_s_*L_ii_)
        K_rt_ = np.sum(L_rt_*L_i_ - L_r_*L_t_*L_ii_)
        K_ss_ = np.sum(L_ss_*L_i_ - L_s_**2*L_ii_)
        K_st_ = np.sum(L_st_*L_i_ - L_s_*L_t_*L_ii_)
        K_tt_ = np.sum(L_tt_*L_i_ - L_t_**2*L_ii_)

        return np.array([[K_rr_, K_rs_, K_rt_],
                         [K_rs_, K_ss_, K_st_],
                         [K_rt_, K_st_, K_tt_]])

    # Define K_x and its derivatives.

    def L_x(r):
        return np.exp(r)*p_x + q_x

#    def L_xr(r):
#        return np.exp(r)*p_x

#    def L_xrr(r):
#        return np.exp(r)*p_x

    def K_x(r):
        return np.sum(np.log(L_x(r)))

#    def dK_x(r):
#        return np.sum(L_xr(r)/L_x(r))

#    def d2K_x(r):
#        return np.sum(L_xrr(r)/L_x(r) - (L_xr(r)**2)/(L_x(r)**2))

    def dK_x_(r):

        a = np.exp(r)*p_x
        b = a + q_x
        c = a/b

        return np.sum(c)

    def d2K_x_(r):

        a = np.exp(r)*p_x
        b = a + q_x
        c = a/b

        return np.sum(c-c**2)

    # Define K_y and its derivatives.

    def L_y(s):
        return np.exp(s)*p_y + q_y

#    def L_ys(s):
#        return np.exp(s)*p_y

#    def L_yss(s):
#        return np.exp(s)*p_y

    def K_y(s):
        return np.sum(np.log(L_y(s)))

#    def dK_y(s):
#        return np.sum(L_ys(s)/L_y(s))

#    def d2K_y(s):
#        return np.sum(L_yss(s)/L_y(s) - (L_ys(s)**2)/(L_y(s)**2))

    def K_y_(s):
        return np.sum(np.log(np.exp(s)*p_y + q_y))

    def dK_y_(s):

        a = np.exp(s)*p_y
        b = a + q_y
        c = a/b

        return np.sum(c)

    def d2K_y_(s):

        a = np.exp(s)*p_y
        b = a + q_y
        c = a/b

        return np.sum(c-c**2)

    # Solve K_xr(r_h) = m_x.

    def dK_x_wrapper(r):
        return np.array([dK_x_(r)-m_x])

    def d2K_x_wrapper(r):
        return np.array([d2K_x_(r)])

    r_h, = fsolve(dK_x_wrapper, 1.0, fprime=d2K_x_wrapper)

    # Solve K_ys(s_h) = m_y.

    def dK_y_wrapper(s):
        return np.array([dK_y_(s)-m_y])

    def d2K_y_wrapper(s):
        return np.array([d2K_y_(s)])

    s_h, = fsolve(dK_y_wrapper, 1.0, fprime=d2K_y_wrapper)

    # Solve dK(r_t, s_t, t_t) = (m_x, m_y, z-0.5).

    def dK_wrapper((r, s, t)):
        return dK_(r, s, t)-m

    def d2K_wrapper((r, s, t)):
        return d2K_(r, s, t)

    r_t, s_t, t_t = fsolve(dK_wrapper, np.ones(3), fprime=d2K_wrapper)

    # Compute saddlepoint approximation.

    w_t = np.sign(t_t)*np.sqrt(2.0*((K_x(r_h)+K_y(s_h))-(m_x*r_h+m_y*s_h) - (K(r_t, s_t, t_t)-(m_x*r_t+m_y*s_t+(z-0.5)*t_t))))
    u_t = 2.0*np.sinh(0.5*t_t)*np.sqrt(np.linalg.det(d2K_(r_t, s_t, t_t))/(d2K_x_(r_h)*d2K_y_(s_h)))

    return norm.cdf(-w_t)-norm.pdf(w_t)*(1.0/w_t-1.0/u_t)
