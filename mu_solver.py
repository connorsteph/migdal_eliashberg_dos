# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:16:36 2017

@author: Connor
"""
import numpy as np
from numpy import sqrt
import tc_func as tf
from zeta_solver import zeta_solver
from chi_solver import chi_solver
import matplotlib.pyplot as plt
import n_sum
from scipy.optimize import brentq
emax = tf.e_max
emin = tf.e_min
epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos
"""
given tc, lambda, w_e, and n, calculates a self consistent mu,
chi, and zeta.
"""

def mu_root_eqn(mu, t, llam, w_e, n, init_chi,
                init_zeta, Nc, maxiter=10,tol=1e-3, damp=0.9, iprint=False):
    dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = sqrt(llam*w_e/2/dos_mu)
    new_zeta = zeta_solver(t, g, w_e, mu,  dos_mu, init_chi,
                           init_zeta, Nc, iprint=False)
    new_chi = chi_solver(t, g, w_e, mu, dos_mu, init_chi,
                         new_zeta, Nc, damp = damp, iprint=False)
    for i in range(1,maxiter+1):
        old_zeta = new_zeta
        old_chi = new_chi
        new_zeta = zeta_solver(t, g, w_e, mu, dos_mu, old_chi, old_zeta,  Nc, iprint=False)
        new_chi = chi_solver(t, g, w_e, mu, dos_mu, old_chi, new_zeta, Nc, damp = damp,
                             iprint=False)
        if (tf.f_compare(new_chi, old_chi) <= tol and
            tf.f_compare(new_chi, old_chi) <= tol):
            break
        if i == maxiter:
            if iprint:
                print('zeta and chi did not converge for mu at t = %g' % t)
                print('damp = %g' % damp)
                print('last difference: %g' % tf.f_compare(new_chi, old_chi))
    return n - n_sum.n_occ(t, g, w_e, tf.dee, emin, emax, mu, dos,
                       new_zeta, new_chi)


def mu_solver(t, llam, w_e, n, init_chi, init_zeta, Nc, tol=1e-3,
              maxiter=5, damp = 0.9, iprint=False):
    dom_lim = len(init_chi)
#    num = 10
#    mu_domain = np.linspace(-0.01, 0.064, num)
#    y = np.zeros(num)
#    for c, mu in enumerate(mu_domain, 0):
#        y[c] = (mu_root_eqn(mu, t, llam, w_e, n, init_chi, init_zeta, Nc, damp=damp))
#    plt.figure()
#    plt.plot(mu_domain, y, 'o-')
#    plt.xlabel('mu')
##    plt.legend(loc = 'best')
#    plt.show()
    mu = brentq(mu_root_eqn, -0.01, 0.01, args=(
            t, llam, w_e, n, init_chi, init_zeta, Nc), xtol=1e-5)
#    print('mu is %g' % mu)
    dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = sqrt(llam*w_e/2/dos_mu)
    new_zeta = zeta_solver(t, g, w_e, mu, dos_mu, init_chi,
                           init_zeta, dom_lim, iprint=False)
    new_chi = chi_solver(t, g, w_e, mu, dos_mu, init_chi,
                         new_zeta, dom_lim, iprint=False)
    for i in range(1,maxiter+1):
        old_zeta = new_zeta
        old_chi = new_chi
        new_zeta = zeta_solver(t, g, w_e, mu, dos_mu, old_chi,
                               old_zeta, dom_lim)
        new_chi = chi_solver(t, g, w_e, mu, dos_mu, old_chi,
                             new_zeta, dom_lim)
        if (tf.f_compare(new_chi, old_chi) <= tol
            and tf.f_compare(new_chi, old_chi) <= tol):
            break
        if i == maxiter:
            if iprint:
                print('zeta and chi not converged for mu')
    return mu, new_zeta, new_chi