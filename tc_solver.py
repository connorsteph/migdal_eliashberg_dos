# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:44:05 2017

@author: Connor
"""

import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq
import tc_func as tf
from mu_solver import mu_solver
import phi_matrix

emax = tf.e_max
emin = tf.e_min
epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos
def tc_root_eqn(t, llam, w_e, n, dom_lim):
    Nc = dom_lim + 20
    init_chi = [tf.init_chi(w) for w in tf.freq_array(1, Nc, t)]
    init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, Nc, t)]
    mu, zeta, chi = mu_solver(t, llam, w_e, n, init_chi, init_zeta, dom_lim)
    dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = np.sqrt(llam*w_e/2/dos_mu)
    p_matrix = phi_matrix.phi_matrix(t, g, w_e, mu, tf.dee, emin, emax,
                                     dos_mu, zeta, chi, dos)
    return np.linalg.det(p_matrix-np.identity(Nc))

def tc_solver(llam, w_e, n, dom_lim, maxiter=100, p_damp=0.3,
               iprint=False, tol=1e-5, damp=0.3, p_tol=1e-2,
               t_tol=5e-2, tc=None,
               ):
    Nc = dom_lim + 30
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
            llam, w_e, n, dom_lim))
    init_chi = [tf.init_chi(w) for w in tf.freq_array(1, Nc, tc)]
    init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, Nc, tc)]
    mu, zeta, chi = mu_solver(tc, llam, w_e, n, init_chi, init_zeta, dom_lim)
    dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = np.sqrt(llam*w_e/2/dos_mu)
    p_matrix = phi_matrix.phi_matrix(tc, g, w_e, tf.dee, emin, emax, dos_mu, zeta, dos)
    eigvals, eigvects = eig(p_matrix)
    roots = [abs(i-1) for i in eigvals]
    phi_v = eigvects[:, np.argmin(roots)]
    phi = np.copy([p/phi_v[0] for p in phi_v])
#    phi = tf.interpolater(tf.freq_array(1, dom_lim, tc), phi_v[:dom_lim])
    return tc, mu, chi, zeta, phi
