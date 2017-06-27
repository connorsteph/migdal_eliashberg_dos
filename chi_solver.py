# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:39:07 2017

@author: Connor
"""

import numpy as np
import matplotlib.pyplot as plt
import tc_func as tf
import chi_sum
emax = tf.e_max
emin = tf.e_min

epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos


def chi_solver(t, g, w_e, mu, dos_mu, dom_lim, zeta_v, init_chi_v, maxiter=150,
                tol=1e-3, iprint=False, damp=0.3):
    Nc = dom_lim + 30
    if iprint:
        plt.figure()
        plt.grid(True)
        plt.plot(tf.m_array(1, dom_lim),
                 [1/w**2 for w in tf.freq_array(1, dom_lim, t)],
                  '.', markersize='2')
    diff_vec = np.empty(maxiter+1)
    new_chi = np.zeros(Nc)
    new_chi = chi_sum.chi(t, g, w_e, tf.dee, emin, emax, mu,
                                  dos_mu, dos, zeta_v, init_chi_v)
    if iprint:
        plt.plot(tf.m_array(1, dom_lim),
                 new_chi[:dom_lim], '--', label='it 0')

    diff_vec[0] = (tf.f_compare(tf.freq_array(1, Nc, t), new_chi))

    """
    we now have a zeta function for the RHS, so we continue to iterate
    until the difference becomes small
    """

    for i in range(1, maxiter+1):
        old_chi = new_chi
        new_chi = chi_sum.chi(
                t, g, w_e, tf.dee, emin, emax, mu, dos_mu, damp,
                dos, zeta_v, old_chi, Nc, tf.nee)
        diff_vec[i] = (tf.f_compare(old_chi, new_chi))

        if(np.mod(i, maxiter // 5) == 0):
            if iprint:
                print('Difference in iteration %i is ' % i, diff_vec[i])
                plt.plot(tf.m_array(1, dom_lim), new_chi[:dom_lim], '-.',
                         label='it %i' % i)
        if (diff_vec[i] < tol):
            if iprint:
                print('chi converged to tol in %i iterations' % i)
                plt.plot(tf.m_array(1, dom_lim), new_chi[:dom_lim], '-.',
                         label='it %i' % i)
                diff_vec = diff_vec[:i+1]
            break
        if i == maxiter:
            print('chi did not converge in %i iterations' % i,
                  'damp = %2.1f' % damp)
            print('last difference: ', diff_vec[-1])

    if iprint:
        print('last difference: ', diff_vec[-1])
        plt.legend(loc='best')
        plt.title('Chi fnc. Damping = %2.2f' % damp)
        plt.savefig('chi_func.pdf', bbox_inches='tight')
        plt.ylabel(r'$\frac{\Chi}{m}$')
        plt.xlabel('m')
        plt.show()

        plt.figure()
        plt.plot(np.log(diff_vec), 'o', markersize=2)
        plt.xlabel('Iteration n')
        plt.ylabel('Log Diff')
        plt.title('Log difference in iterated fnc. Damping =%2.2f' % damp)
        plt.show()

    return new_chi[:dom_lim]
