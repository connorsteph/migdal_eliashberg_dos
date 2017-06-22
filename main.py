# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:40:00 2017

@author: Connor
"""

import numpy as np
from tc_calc import tc_calc
import tc_func as tf
import phi_sum
#print(phi_sum.phi_sum_init.__doc__)

emin = tf.e_min
emax = tf.e_max
"""
Problem params
"""
D = 1/(emax-emin)
w_e = 16/2.5
lam_want = 2
g = np.sqrt(lam_want*w_e/2/tf.dos(0))

"""
Algorithm params
"""

dom_lim = 200
p_damp = 0.3
maxiter = 50
tol = 1e-5
damp = 0.3

#domain = np.linspace(0.5, 20, 10)
#plt.figure()
#plt.grid(True)
#plt.plot(domain, [tc_calc(np.sqrt(lam_want*(Q*w_e)*w_e/2), w_e, Q*w_e,
#                          dom_lim, maxiter=maxiter,
#                          tol=tol, p_tol=5e-5, t_tol=5e-2, plot=False)/w_e
#                  for Q in domain])
#plt.title(r'$T_c$ as a function of $\frac{D}{\omega_E}$', fontsize=22)
#plt.ylabel(r'$\frac{T_c}{\omega_E}$', fontsize=18)
#plt.xlabel(r'$\frac{D}{\omega_E}$', fontsize=18)
#plt.savefig('tc_plot_lam_%g.pdf' % lam_want, bbox_inches='tight')


tc_calc(g, w_e, D,
        dom_lim, maxiter=maxiter,
        tol=tol, p_tol=5e-5, t_tol=5e-2, plot=False, iprint=False)



#num = 100
#llam = 2*g**2/(D*w_e)
#zeta = zeta_solver(t, g, w_e, num, D, damp=0.3, maxiter=150, iprint=False)
#plt.figure()
#plt.grid(True)
#domain = tf.m_array(1, num)
#for m in [20, 50, 90]:
#    y = [init_summand(w_n, tf.freq_m(m, t), init_phi, zeta, w_e, t, D)
#         for w_n in tf.freq_array(1, num, t)]
#    plt.plot(domain, y, label=('m = %i' % m))
#plt.legend(loc='best')
#plt.xlabel('n')
##
#plt.figure()
#plt.grid(True)
#domain = np.arange(1, 200, 10)
#for m in [1]:
#    y = [llam*np.pi*t*(tf.matsu_sum(1, M, t, init_summand, tf.freq_m(m, t),
#                       init_phi, zeta, w_e, t, D)) - init_phi(m, t)
#         for M in domain]
#    plt.plot(domain, y, label=('m = %g' % m))
#plt.legend(loc='best')


#plt.figure()
#domain = np.linspace(1, 21)
#for q in [100*s, 200*s]:
#    y=[]
#    for t in domain:
#        zeta = zeta_solver(t, np.sqrt(q), w_e, 50, D,
#                                     damp=0.3, maxiter=150, iprint=False)
#        y.append(zeta(tf.freq_m(1, t))/tf.freq_m(1, t))
#    plt.plot(domain, y, label=('q = %g' % (q)))
#plt.plot(domain, [(2*(200*s)/(D*w_e))+1 for t in domain], label='correct')
##plt.ylim([0, 1.5])
#plt.xlabel('t/w_e')
#plt.legend(loc='best')
#plt.grid(True)
#zeta = zeta_solver(t, g, w_e, 150, D, damp=0.3, maxiter=150, iprint=True)
