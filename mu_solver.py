# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:16:36 2017

@author: Connor
"""
import numpy as np
from numpy import sqrt
import tc_func as tf
#from zeta_solver import zeta_solver
#from chi_solver import chi_solver
import matplotlib.pyplot as plt
import n_sum
from scipy.optimize import brentq
import mu_zeta_chi
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
                init_zeta, Nc, damp, maxiter,tol=1e-3, iprint=False):
#    print(n)
#    print(t/w_e)
    if np.sign(mu) < 0:
        for c, val in enumerate(init_chi):
            init_chi[c] = -val
    dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = sqrt(llam*w_e/2/dos_mu)
    zeta, chi, _ = mu_zeta_chi.simult_mu_func(t,g,w_e,mu,dos_mu,emin,emax,tf.dee,damp,
                                           maxiter,tol,dos,init_zeta,init_chi)
#    print(n)
    return n - n_sum.n_occ(t, g, w_e, tf.dee, emin, emax, mu, dos,
                       zeta, chi)
#    return mu_zeta_chi.mu_root_eqn(t,g,w_e,tf.dee,emin,emax,mu,dos_mu,n,maxiter,
#                                   damp,tol,dos,init_zeta,init_chi)


#    new_chi = chi_solver(
#            t, g, w_e, mu, dos_mu, init_chi, init_zeta, Nc,
#            maxiter=100, damp = damp, iprint=True)
#    new_zeta = zeta_solver(t, g, w_e, mu,  dos_mu, new_chi,
#                           init_zeta, Nc, iprint=False)
#
#    for i in range(1,maxiter+1):
#        old_zeta = new_zeta
#        old_chi = new_chi
#        new_zeta = zeta_solver(t, g, w_e, mu, dos_mu, old_chi,
#                               old_zeta,  Nc, iprint=False)
#        new_chi = chi_solver(t, g, w_e, mu, dos_mu, old_chi,
#                             new_zeta, Nc, damp = damp,
#                             iprint=False, maxiter=100)
#
#        if (tf.f_compare(new_chi, old_chi) <= tol and
#            tf.f_compare(new_chi, old_chi) <= tol):
#            break
#        if i == maxiter:
#            if iprint:
#                print('zeta and chi did not converge for mu at t = %g' % t)
#                print('damp = %g' % damp)
#                print('last difference: %g' % tf.f_compare(new_chi, old_chi))
#    return n - n_sum.n_occ(t, g, w_e, tf.dee, emin, emax, mu, dos,
#                       new_zeta, new_chi)


def mu_solver(t, llam, w_e, n, init_chi, init_zeta, Nc, tol=1e-3,
              maxiter=150, damp = 0.9, iprint=False):
    cutoff = 0.0
#    print('flag1')
#    num = 10
#    mu_domain = np.linspace(emin+cutoff, emax-cutoff, num)
#    y = np.zeros(num)
#    for c, mu in enumerate(mu_domain, 0):
#        y[c] = (mu_root_eqn(mu, t, llam, w_e, n, init_chi,
#         init_zeta, Nc, damp, maxiter))
#    plt.figure()
#    plt.plot(mu_domain, y, 'o-')
#    plt.xlabel('mu')
#    plt.show()
    mu = brentq(mu_root_eqn, emin+cutoff, emax-cutoff, args=(
            t, llam, w_e, n, init_chi, init_zeta, Nc, damp, maxiter), xtol=1e-5)
#    print('flag2')
    dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = sqrt(llam*w_e/2/dos_mu)
    zeta, chi, _ = mu_zeta_chi.simult_mu_func(t,g,w_e,mu,dos_mu,emin,emax,tf.dee,damp,
                                           maxiter,tol,dos,init_zeta,init_chi)
#    print(t/w_e)


    #print('mu is %g' % mu)
    return mu, zeta, chi