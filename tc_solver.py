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
from omp_phi_matrix import phi_mat

emax = tf.e_max
emin = tf.e_min
dos = tf.dos
dos_avg = tf.dos_avg

def tc_root_eqn(t,  llam, w_e, n, Nc, maxiter, damp, tol):
    mu, zeta, chi, z = mu_solver(t, llam, w_e, n, Nc, maxiter=maxiter,
                                 damp=damp, tol=1e-8)
    print('First index of chi:', chi[0])
    p_matrix = phi_mat.phi_matrix(t, llam, w_e, mu, dos_avg, tf.dee, emin,
                                      zeta, chi, dos, tf.q_list, tf.p_list)
    print(p_matrix)
    val = np.linalg.det(p_matrix - np.identity(Nc))
    print(val, t/w_e)
    return val

def tc_solver(llam, w_e, n, dom_lim, maxiter=150,
              mu_tol=1e-5, damp=0.2):
    print("Solving for tc")
    Nc = dom_lim
    min_bound = 0.1 * w_e
    max_bound = 0.3 * w_e
    tc = brentq(tc_root_eqn, min_bound, max_bound, args=(
        llam, w_e, n, Nc, maxiter, damp, mu_tol), xtol=1e-10)
    mu, zeta, chi, z = mu_solver(tc, llam, w_e, n, Nc, maxiter=maxiter, damp=damp, tol=1e-8)
    p_matrix = phi_mat.phi_matrix(tc, llam, w_e, mu, dos_avg, tf.dee, emin,
                                     zeta, chi, dos, tf.q_list, tf.p_list)
    eigvals, eigvects = eig(p_matrix)
    roots = [abs(i - 1) for i in eigvals]
    phi_v = eigvects[:, np.argmin(roots)]
    phi = np.copy([p / phi_v[0] for p in phi_v])
    return [tc, mu, chi[:dom_lim], zeta[:dom_lim], phi[:dom_lim], z[:dom_lim]]
