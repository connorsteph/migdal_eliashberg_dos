# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:44:05 2017

@author: Connor
"""

import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
from scipy.optimize import brentq
import tc_func as tf
from zeta_solver import zeta_solver
import phi_sum

emax = tf.e_max
emin = tf.e_min
epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos


def tc_root_eqn(
        t, g, w_e, mu, dos_mu, phi, dom_lim):
    Nc = 25
    zeta_v = zeta_solver(t, g, w_e, mu, dos_mu, Nc)
    llam = 2*dos_mu*g**2/w_e
    phi_v = [phi(w) for w in tf.freq_array(1, Nc, t)]
    return np.pi*llam/dos_mu*t*phi_sum.phi_sum_init(
                t, g, w_e, tf.freq_m(1, t), tf.dee, emin, emax, phi_v, zeta_v,
                dos, Nc, tf.nee) - phi_v[0]

def phi_solver(g, w_e, mu, dos_mu, dom_lim, init_phi, maxiter=100, p_damp=0.3,
               iprint=False, tol=1e-5, damp=0.3, p_tol=1e-2,
               t_tol=5e-2, tc=None,
               ):
    l_root = 0.01*w_e
    diff_vec = []
    Nc = dom_lim + 25
    new_phi = np.zeros(Nc)
    if tc is None:
        if iprint:
            print('Solving for initial tc')
            plt.figure()
            plt.grid(True)
            num = 11
            t_domain = np.linspace(l_root, w_e, num)
            y = np.zeros(num)
            for c, t in enumerate(t_domain, 0):
                y[c] = (tc_root_eqn(t, g, w_e, mu, dos_mu, init_phi, dom_lim))
            plt.plot(t_domain, y, 'o-')
            plt.xlim([0, w_e])
            plt.xlabel('t')
            plt.show()
        tc = brentq(tc_root_eqn, l_root, w_e, args=(
            g, w_e, mu, dos_mu, init_phi, dom_lim))
    if iprint:
        print('Tc/w_e = %5.4g' % (tc/w_e))
    zeta_v = zeta_solver(tc, g, w_e, mu, dos_mu, Nc, tol=tol, iprint=False)
    init_phi_v = [init_phi(w) for w in tf.freq_array(1, Nc, tc)]
    if iprint:
        plt.figure()
        plt.grid(True)
        plt.ylim([0, 1])
        plt.xlabel('w_m')
        plt.plot([w/w_e for w in tf.freq_array(1, dom_lim, tc)],
                 [init_phi_v[i]/init_phi_v[0]
                 for i in tf.m_array(1, dom_lim)], label='initial')

    if iprint:
        print('iterating initial phi')
    new_phi = phi_sum.phi_init( tc, g, w_e, tf.dee, emin, emax, dos_mu,
                               init_phi_v, zeta_v, dos, Nc, tf.nee)

    if iprint:
        print('converging phi')
    for i in range(1, maxiter+1):
        old_phi = new_phi
        new_phi = phi_sum.phi(tc, g, w_e, tf.dee, emin, emax, dos_mu,
                              p_damp, old_phi, zeta_v, dos, Nc, tf.nee)
        diff_vec.append(tf.f_compare(old_phi, new_phi))

        if iprint:
            if(np.mod(i, maxiter // 5) == 0):
                print('Difference in iteration %i is ' % i, diff_vec[i-1])
                print('pdamp = %g' % p_damp)
                plt.plot([w/w_e for w in tf.freq_array(1, dom_lim, tc)],
                         new_phi[:dom_lim], '-.', label='it %i' % i)
        if (diff_vec[i-1] < p_tol):
            if iprint:
                print('phi converged to p_tol in %i iterations' % i)
                print('p_tol = %3.2g, Nc = %i, g = %3.2g, w_e = %3.2g,\n\
                      p_damp = %3.2g'
                      % (p_tol, Nc, g, w_e, p_damp))
                plt.plot([w/w_e for w in tf.freq_array(1, dom_lim, tc)],
                         new_phi[:dom_lim], '-.', label='it %i' % i)
            break

    if iprint:
        plt.legend(loc='best')
        plt.savefig('phi_iter.pdf', bbox_inches='tight')
        plt.show()

        plt.figure()
        plt.plot(np.log(diff_vec), 'o', markersize=2)
        plt.xlabel('Iteration n')
        plt.ylabel('Log Diff')
        plt.title('Log difference in phi iter. Damping = %2.2f' % p_damp)
        plt.show()

    phi = tf.interpolater(tf.freq_array(1, dom_lim, tc), new_phi[:dom_lim])
    try:
        old_vals = [init_phi(w)/init_phi(tf.freq_m(1, tc))
                    for w in tf.freq_array(1, dom_lim, tc)]
    except TypeError:
        old_vals = [init_phi(tf.matsu_index(w, tc), tc)/init_phi(1, tc)
                    for w in tf.freq_array(1, dom_lim, tc)]

    diff = tf.f_compare(old_vals, new_phi[:dom_lim])
    new_tc = brentq(tc_root_eqn, l_root, w_e, args=(
            g, w_e, mu, dos_mu, phi, dom_lim))
    t_diff = abs(tc - new_tc)

    dom_lim = np.int(tc/new_tc*(dom_lim-1/2)+1/2)
    """
    keeps new dom_lim from extending past domain of the initial phi,
    means that the initial dom_lim needs to be larger to compensate
    if tc increases
    """
    if (diff < p_tol and t_diff < t_tol):
        return phi, new_tc, dom_lim
    else:
        return phi_solver(g, w_e, mu, dos_mu, dom_lim, phi,
                          maxiter=maxiter, p_damp=p_damp,
                          iprint=iprint, tol=tol, tc=new_tc
                          )
