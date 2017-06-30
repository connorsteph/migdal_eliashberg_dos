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
import matplotlib.pyplot as plt

emax = tf.e_max
emin = tf.e_min
epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos
dos_avg = tf.dos_avg

def tc_root_eqn(t, llam, w_e, n, Nc, maxiter, damp=0.6):

#    t_samp.append(t)
#    print(t/w_e)
    init_chi = [tf.init_chi(w) for w in tf.freq_array(1, Nc, t)]
    init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, Nc, t)]
#    print('flag0')
    mu, zeta, chi, z = mu_solver(t, llam, w_e, n, init_chi, init_zeta,
                              Nc, maxiter=maxiter, damp=damp)
#    plt.figure()
#    plt.plot(chi)
#    plt.show()
#    print(mu)
#    print(t/w_e)
    p_matrix = phi_matrix.phi_matrix(t, llam, w_e, mu, dos_avg, tf.dee, emin, emax,
                                     zeta, chi, dos)
    return np.linalg.det(p_matrix-np.identity(Nc))

def tc_solver(llam, w_e, n, dom_lim, maxiter=150, p_damp=0.3,
               iprint=False, tol=1e-5, damp=0.6, p_tol=1e-2, tc=None,
               ):
#    print('tc solver')
    Nc = dom_lim + 20
    min_bound = 0.1*w_e
    max_bound = 0.9*w_e
#    num = 10
#    t_domain = np.linspace(min_bound, max_bound, num)
#    y = np.zeros(num)
#    for c, t in enumerate(t_domain, 0):
#        y[c] = (tc_root_eqn(t, llam, w_e, n, dom_lim, maxiter))
#    plt.figure()
#    plt.grid(True)
#    plt.plot([t/w_e for t in t_domain], y, 'o-',label = 'w_e = %5.4g' % w_e)
#    plt.xlabel('t')
#    plt.legend(loc = 'best')
#    plt.show()
    tc = brentq(tc_root_eqn, min_bound, max_bound, args=(
            llam, w_e, n, Nc, maxiter), xtol=1e-2)
    init_chi = [tf.init_chi(w) for w in tf.freq_array(1, Nc, tc)]
    init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, Nc, tc)]
    mu, zeta, chi, z = mu_solver(tc, llam, w_e, n, init_chi, init_zeta, Nc)
    p_matrix = phi_matrix.phi_matrix(tc, llam, w_e, mu, dos_avg, tf.dee, emin, emax,
                                      zeta, chi, dos)
    eigvals, eigvects = eig(p_matrix)
    roots = [abs(i-1) for i in eigvals]
    phi_v = eigvects[:, np.argmin(roots)]
    phi = np.copy([p/phi_v[0] for p in phi_v])
    return [tc, mu, chi[:], zeta[:], phi[:], z[:]]
