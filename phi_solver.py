# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:44:05 2017

@author: Connor
"""

import numpy as np
from scipy.linalg import eig
#from scipy.integrate import quad
from matplotlib import pyplot as plt
from scipy.optimize import brentq
import tc_func as tf
from zeta_solver import zeta_solver
import phi_matrix
import phi_sum

emax = tf.e_max
emin = tf.e_min
epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos
def tc_root_eqn(t, g, w_e, mu, dos_mu, dom_lim):
    zeta_v = zeta_solver(t, g, w_e, mu, dos_mu, dom_lim+20)
    p_matrix = phi_matrix.phi_matrix(t, g, w_e, tf.dee, emin, emax,
                                     dos_mu, zeta_v, dos)
    return np.linalg.det(p_matrix-np.identity(dom_lim+20))

def phi_solver(g, w_e, mu, dos_mu, dom_lim, init_phi, maxiter=100, p_damp=0.3,
               iprint=False, tol=1e-5, damp=0.3, p_tol=1e-2,
               t_tol=5e-2, tc=None,
               ):
    l_root = 0.01*w_e
#    plt.figure(0)
#    plt.clf
#    plt.grid(True)
#    num = 10
#    t_domain = np.linspace(l_root, w_e, num)
#    y = np.zeros(num)
#    for c, t in enumerate(t_domain, 0):
#        y[c] = (tc_root_eqn(t, g, w_e, mu, dos_mu, dom_lim))
#    plt.plot([t/w_e for t in t_domain], y, 'o-',label = 'w_e = %5.4g' % w_e)
#    plt.xlim([0, 1])
#    plt.xlabel('t')
#    plt.legend(loc = 'best')
#    plt.show()
    tc = brentq(tc_root_eqn, l_root, w_e, args=(
            g, w_e, mu, dos_mu, dom_lim))

    zeta_v = zeta_solver(tc, g, w_e, mu, dos_mu, dom_lim)
    p_matrix = phi_matrix.phi_matrix(tc, g, w_e, tf.dee, emin, emax, dos_mu, zeta_v, dos)
    eigvals, eigvects = eig(p_matrix)
    roots = [abs(i-1) for i in eigvals]
    phi_v = eigvects[:, np.argmin(roots)]
    phi_v = np.copy([p/phi_v[0] for p in phi_v])
    phi = tf.interpolater(tf.freq_array(1, dom_lim, tc), phi_v[:dom_lim])
    return tc, phi, phi_v
