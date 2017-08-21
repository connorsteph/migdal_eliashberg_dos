# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:16:36 2017

@author: Connor
"""
import numpy as np
import tc_func as tf
from scipy.optimize import brentq
import omp_mu_zeta_chi
import matplotlib.pyplot as plt
emax = tf.e_max
emin = tf.e_min
dos = tf.dos
dos_avg = tf.dos_avg

"""
given tc, lambda, w_e, and n, calculates a self consistent mu,
chi, and zeta.
"""

def mu_root_eqn(mu, t, llam, w_e, n, Nc, damp, maxiter, tol):
    # if np.sign(mu) < 0:
    #     for c, val in enumerate(init_chi):
    #         init_chi[c] = -val
    val, iterations = omp_mu_zeta_chi.mu_equations.mu_root_eqn(
        t, llam, w_e, mu, dos_avg, n, emin, tf.dee, damp, maxiter,
        tol, dos, tf.q_list, tf.p_list, Nc, tf.nee)
    if iterations < 0:
        print('slow convergence, damp = %g, mu = %g, t/w_e = %g' %
              (damp, mu, t / w_e))
    print("%14.8g  | %7.5g" %(val, mu))
    return val


def mu_solver(t, llam, w_e, n, Nc, tol=1e-18,
              maxiter=150, damp=0.2):
    print("Solving for mu with tc/we %g" %(t/w_e))
    print("        Val            Mu ")   
    cutoff = 0.0
    mu = brentq(mu_root_eqn, emin + cutoff, emax - cutoff, args=(
        t, llam, w_e, n, Nc, damp, maxiter, tol), xtol=1e-5)
    zeta, chi, _ = omp_mu_zeta_chi.mu_equations.simult_mu_func(
        t, llam, w_e, mu, dos_avg, emin, tf.dee, damp,
        maxiter, tol, dos, tf.q_list, tf.p_list, Nc)
    z = [zeta[i] / tf.freq_m(i+1, t) for i in range(len(zeta))]
    return mu, zeta, chi, z
