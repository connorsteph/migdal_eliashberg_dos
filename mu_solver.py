# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:16:36 2017

@author: Connor
"""
from numpy import sqrt
import tc_func as tf
from zeta_solver import zeta_solver
from chi_solver import chi_solver
import n_sum
from scipy.optimize import brentq
emax = tf.e_max
emin = tf.e_min
epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos
"""
given tc, lambda, w_e, and n, calculates a self consistent mu, chi, and zeta.
"""

def mu_root_eqn(mu, t, llam, w_e, n, init_chi_v, init_zeta_v,
                     dos_mu, dom_lim, maxiter=50,tol=1e-3, damp=0.3):
    dos_mu = dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = sqrt(llam*w_e/2/dos_mu)
    new_zeta = zeta_solver(t, g, w_e, mu, init_chi_v, init_zeta_v, dos_mu, dom_lim+20)
    new_chi = chi_solver(t, g, w_e, mu, init_chi_v, new_zeta, dos_mu, dom_lim+20)
    for i in range(1,maxiter+1):
        old_zeta = new_zeta
        old_chi = new_chi
        new_zeta = zeta_solver(t, g, w_e, mu, old_chi, old_zeta, dos_mu, dom_lim+20)
        new_chi = chi_solver(t, g, w_e, mu, old_chi, new_zeta, dos_mu, dom_lim+20)
        if (tf.f_compare(new_chi, old_chi) <= tol and tf.f_compare(new_chi, old_chi) <= tol):
            break
        if i == maxiter:
            print('zeta and chi not converged for mu')
    return n - n_sum.n_occ(t, g, w_e, emin, emax, mu, dos,
                       new_zeta, new_chi)


def mu_solver(t, llam, w_e, n, init_chi, init_zeta,
              dos_mu, tol=1e-3, maxiter=50):
    dom_lim = len(init_chi)
    mu = brentq(mu_root_eqn, emin, emax, args=(
            t, llam, w_e, n, init_chi, init_zeta, dos_mu, dom_lim))
    dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
    g = sqrt(llam*w_e/2/dos_mu)
    new_zeta = zeta_solver(t, g, w_e, mu, init_chi, init_zeta, dos_mu, dom_lim+20)
    new_chi = chi_solver(t, g, w_e, mu, init_chi, new_zeta, dos_mu, dom_lim+20)
    for i in range(1,maxiter+1):
        old_zeta = new_zeta
        old_chi = new_chi
        new_zeta = zeta_solver(t, g, w_e, mu, old_chi, old_zeta, dos_mu, dom_lim+20)
        new_chi = chi_solver(t, g, w_e, mu, old_chi, new_zeta, dos_mu, dom_lim+20)
        if (tf.f_compare(new_chi, old_chi) <= tol and tf.f_compare(new_chi, old_chi) <= tol):
            break
        if i == maxiter:
            print('zeta and chi not converged for mu')
    return mu, new_zeta, new_chi