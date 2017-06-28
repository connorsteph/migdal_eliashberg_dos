# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:40:00 2017

@author: Connor
"""

import numpy as np
from tc_calc import tc_calc
import tc_func as tf
from time import time
import matplotlib.pyplot as plt
from tc_solver import tc_solver
start = time()
emin = tf.e_min
emax = tf.e_max
"""
Problem params
"""
ttp = tf.ttp
D = (emax-emin)
w_e = 16/2.5*tf.ttp
lam_want = 2*tf.ttp
n = 1.0
print('Lambda = %g' %lam_want)
"""
Algorithm params
"""

dom_lim = 25
p_damp = 0.3
maxiter = 50
tol = 1e-5
damp = 0.9


tc, mu, chi, zeta, phi = tc_solver(
        lam_want, w_e, n, dom_lim, maxiter=maxiter,
        tol=tol, p_tol=tol, t_tol=5e-2, iprint=False, damp=damp)
print('Tc/we is %g' % (tc/w_e))


#plt.figure()
#plt.plot(tf.dos)
#points = 20
#domain = np.linspace(0.01*ttp, 1.0*ttp, points)
#rrange = np.empty(points)
#for c, w in enumerate(domain):
#    g = np.sqrt(lam_want*w/2/dos_mu)
##    inner_start = time()
##    print((c+1)/points, w/D)
#    rrange[c] = tc_calc(g, w, n, mu, dom_lim)/w
##    print('tc is ', rrange[c]/w)
##    print(time() - inner_start)
#np.savetxt('bcc_dos_range_lam_1.dat', rrange, delimiter=',')
#np.savetxt('bcc_dos_domain_lam_1.dat', domain, delimiter=',')

#bcc_domain_lam_2 = np.loadtxt('bcc_dos_domain_lam_2.dat', delimiter=',')
#bcc_vals_lam_2 = np.loadtxt('bcc_dos_range_lam_2.dat', delimiter=',')
#const_domain_lam_2 = np.loadtxt('const_dos_domain_lam_2.dat', delimiter=',')
#const_vals_lam_2 = np.loadtxt('const_dos_range_lam_2.dat', delimiter=',')
#const_domain_lam_1 = np.loadtxt('const_dos_domain_lam_1.dat', delimiter=',')
#const_vals_lam_1 = np.loadtxt('const_dos_range_lam_1.dat', delimiter=',')
#bcc_domain_lam_1 = np.loadtxt('bcc_dos_domain_lam_1.dat', delimiter=',')
#bcc_vals_lam_1 = np.loadtxt('bcc_dos_range_lam_1.dat', delimiter=',')
#
#plt.figure()
#plt.grid(True)
#plt.plot(const_domain_lam_2, const_vals_lam_2, '-.', label = 'const_dos_lam_2')
#plt.plot(const_domain_lam_1, const_vals_lam_1, '--', label = 'const_dos_lam_1')
#plt.plot(bcc_domain_lam_2, bcc_vals_lam_2, '-.', label = 'bcc_lam_2')
#plt.plot(bcc_domain_lam_1, bcc_vals_lam_1, '--', label = 'bcc_lam_1')
#plt.xlabel(r'$\frac{\omega_E}{t}$', fontsize=14)
#plt.ylabel(r'$\frac{T_C}{\omega_E}$', fontsize=14)
#plt.legend(loc='best')
#plt.title('Critical Temp. versus Char. Freq. ')
#plt.savefig('tc_vs_w_e_const_and_bcc.pdf', bbox_inches='tight', dpi = 150)

#plt.plot([x/D for x in domain], rrange)
#plt.xlabel(r'$\frac{\omega_E}{t}$', fontsize=14)
#plt.ylabel(r'$\frac{T_C}{\omega_E}$', fontsize=14)
#plt.title('Critical Temp. versus Char. Freq. for Const. DOS 1/D')
#plt.savefig('tc_vs_w_e_lam_2_const_dos_dee_5d-2.pdf', bbox_inches='tight', dpi = 300)
#plt.figure()
#plt.plot([1,2,3,4,5], [tc_calc(g, w_e, n, mu, dom_lim, maxiter=maxiter, tol=tol,
#             p_tol=tol, t_tol=5e-2, plot=False, iprint=False)/w_e for tol in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]])
print('Runtime = %g' % (time() - start))
