# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:16:42 2017

@author: Connor
"""

import numpy as np
import matplotlib.pyplot as plt
import tc_func as tf
import zeta_sum
emax = tf.e_max
emin = tf.e_min

epsrel = 1e-4
epsabs = 1e-4
dos = tf.dos


def zeta_solver(t, g, w_e, mu, dos_mu, chi, init_zeta, dom_lim, maxiter=150,
                tol=1e-3, iprint=False, damp=0.3):
#    print('zeta solver')
    """
    attempts to converge the zeta function, from an initial guess
    where zeta(iw_m)=w_m.
    A damping factor has been introduced to aid in convergence

    Input: Tc, g, w_e, limits of frequency indices to be calc'd over,
    limits of Matsu sum, max Ncber of iterations, damping factor

    Output: a callable function zeta(x), converged for the given Tc


    index i of new_zeta is equal to new_zeta evaluated at m=i+1
    """
#    Nc = dom_lim + 20
    if iprint:
        plt.figure()
        plt.grid(True)
        plt.plot(tf.m_array(1, dom_lim), tf.freq_array(
                1, dom_lim, t), '.', markersize='2')
    diff_vec = np.empty(maxiter+1)
    new_zeta = zeta_sum.zeta(t, g, w_e, mu, tf.dee, emin, emax,
                                  dos_mu, damp, dos, init_zeta, chi)
    if iprint:
        plt.plot(tf.m_array(1, dom_lim),
                 new_zeta[:dom_lim], '--', label='it 0')

    diff_vec[0] = (tf.f_compare(tf.freq_array(1, dom_lim, t), new_zeta))

    """
    we now have a zeta function for the RHS, so we continue to iterate
    until the difference becomes small
    """

    for i in range(1, maxiter+1):
        old_zeta = new_zeta
        new_zeta = zeta_sum.zeta(
                t, g, w_e, mu, tf.dee, emin, emax,
                dos_mu, damp, dos, old_zeta, chi)
        diff_vec[i] = (tf.f_compare(old_zeta, new_zeta))

        if(np.mod(i, maxiter // 5) == 0):
            if iprint:
                print('Difference in iteration %i is ' % i, diff_vec[i])
                plt.plot(tf.m_array(1, dom_lim), new_zeta[:dom_lim], '-.',
                         label='it %i' % i)
        if (diff_vec[i] < tol):
            if iprint:
                print('zeta converged to tol in %i iterations' % i)
                plt.plot(tf.m_array(1, dom_lim), new_zeta[:dom_lim], '-.',
                         label='it %i' % i)
                diff_vec = diff_vec[:i+1]
            break
        if i == maxiter:
            print('Zeta did not converge in %i iterations' % i,
                  'damp = %2.1f' % damp)
            print('last difference: ', diff_vec[-1])

    if iprint:
        print('last difference: ', diff_vec[-1])
        plt.legend(loc='best')
        plt.title('Zeta fnc. Damping = %2.2f' % damp)
        plt.savefig('zeta_func.pdf', bbox_inches='tight')
        plt.ylabel('Zeta(m)')
        plt.xlabel('m')
        plt.show()

        plt.figure()
        plt.plot(np.log(diff_vec), 'o', markersize=2)
        plt.xlabel('Iteration n')
        plt.ylabel('Log Diff')
        plt.title('Log difference in iterated fnc. Damping =%2.2f' % damp)
        plt.show()
#    zeta = tf.interpolater(tf.freq_array(1, dom_lim, t), new_zeta[:dom_lim])
#    return zeta, new_zeta[:dom_lim]
    return new_zeta[:dom_lim]
