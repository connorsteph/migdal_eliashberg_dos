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


def chi_solver(t, g, w_e, mu, dos_mu, init_chi, zeta, Nc, maxiter=150,
                tol=1e-3, iprint=False, damp=0.90):
#    print('mu solver')
    if iprint:
        plt.figure(5)
        plt.grid(True)
        plt.plot(tf.m_array(1, Nc),
                init_chi, '.', markersize='2', label='initial')
    diff_vec = np.empty(maxiter+1)
    new_chi = chi_sum.chi(t, g, w_e, mu, tf.dee, emin, emax,
                          dos_mu, damp, dos, zeta, init_chi, Nc, tf.nee)
    if iprint:
        plt.figure(5)
        plt.plot(tf.m_array(1, Nc),
                 new_chi[:Nc], '--', label='it 0')
    diff_vec[0] = (tf.f_compare(tf.freq_array(1, Nc, t), new_chi))

    """
    we now have a zeta function for the RHS, so we continue to iterate
    until the difference becomes small
    """

    for i in range(1, maxiter+1):
        old_chi = new_chi
        new_chi = chi_sum.chi(
                t, g, w_e, mu, tf.dee, emin, emax, dos_mu, damp,
                dos, zeta, old_chi, Nc, tf.nee)
        diff_vec[i] = (tf.f_compare(old_chi, new_chi))
        if(np.mod(i, maxiter // 5) == 0):
            if iprint:
                print('Difference in iteration %i is ' % i, diff_vec[i])
                plt.figure(5)
                plt.plot(tf.m_array(1, Nc), new_chi[:Nc], '-.',
                         label='it %i' % i)
        if (diff_vec[i] < tol):
            if iprint:
                print('chi converged to tol in %i iterations' % i)
                plt.figure(5)
                plt.plot(tf.m_array(1, Nc), new_chi[:Nc], '-.',
                         label='it %i' % i)

                diff_vec = diff_vec[:i+1]
            break
        if i == maxiter:
            print('chi did not converge in %i iterations' % i,
                  'damp = %2.1f, temp = %g' % (damp,t))
            if iprint is False:
                print('last difference: ', diff_vec[-1])

    if iprint:
        plt.figure(5)
        plt.grid(True)
        print('last difference: ', diff_vec[-1])
        plt.legend(loc='best')
        plt.title('Chi fnc. Damping = %2.2f' % damp)
        plt.ylabel(r'$\chi_m$',fontsize=18)
        plt.xlabel('m')
        plt.savefig('chi_func.pdf', bbox_inches='tight')
        plt.show()

        plt.figure()
        plt.plot(np.log(diff_vec), '-o', markersize=2)
        plt.xlabel('Iteration n')
        plt.ylabel('Log Diff')
        plt.title('Log difference in iterated fnc. Damping =%2.2f' % damp)
        plt.show()

    return new_chi[:Nc]
