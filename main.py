# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:40:00 2017

@author: Connor
"""

import numpy as np
from tc_calc import tc_calc
import tc_func as tf
from time import time
from zeta_solver import zeta_solver
import matplotlib.pyplot as plt

start = time()
emin = tf.e_min
emax = tf.e_max
"""
Problem params
"""
D = (emax-emin)
w_e = 16/2.5
lam_want = 2
g = np.sqrt(lam_want*w_e/2/tf.dos[np.int(tf.nee/2 + 1)])
print('Lambda = %g' %lam_want)
"""
Algorithm params
"""

dom_lim = 100
p_damp = 0.3
maxiter = 50
tol = 1e-5
damp = 0.3

points = 10
domain = np.linspace(0.1*D, 1*D, points)
rrange = np.empty(points)
for c, w in enumerate(domain):
    g = np.sqrt(lam_want*w/2/tf.dos[np.int(tf.nee/2 + 1)])
    inner_start = time()
    print((c+1)/points)
    rrange[c] = tc_calc(g, w, D,
        dom_lim, maxiter=maxiter,
        tol=1e-5, p_tol=5e-5, t_tol=5e-2, plot=False, iprint=False)/w
    print(time() - inner_start)
np.savetxt('bcc_dos_range.dat', rrange, delimiter=',')
np.savetxt('bcc_dos_domain.dat', domain, delimiter=',')

#bcc_domain = np.loadtxt('BCC_dos_domain.dat', delimiter=',')
#bcc_vals = np.loadtxt('BCC_dos_range.dat', delimiter=',')
#const_domain = np.loadtxt('const_dos_domain.dat', delimiter=',')
#const_vals = np.loadtxt('const_dos_range.dat', delimiter=',')

#plt.figure()
#plt.grid(True)
#plt.plot([x/D for x in const_domain], const_vals, label = 'const_dos')
#plt.plot([x/D for x in bcc_domain], bcc_vals, label = 'BCC')
#plt.xlabel(r'$\frac{\omega_E}{D}$', fontsize=14)
#plt.ylabel(r'$\frac{T_C}{\omega_E}$', fontsize=14)
#plt.legend(loc='best')
#plt.title('Critical Temp. versus Char. Freq. ')
#plt.savefig('tc_vs_w_e_lam_2_const_dos_dee_5d-2.pdf', bbox_inches='tight', dpi = 300)

#plt.plot([x/D for x in domain], rrange)
#plt.xlabel(r'$\frac{\omega_E}{D}$', fontsize=14)
#plt.ylabel(r'$\frac{T_C}{\omega_E}$', fontsize=14)
#plt.title('Critical Temp. versus Char. Freq. for Const. DOS 1/D')
#plt.savefig('tc_vs_w_e_lam_2_const_dos_dee_5d-2.pdf', bbox_inches='tight', dpi = 300)

#tc = tc_calc(g, w_e, D,
#        dom_lim, maxiter=maxiter,
#        tol=tol, p_tol=5e-5, t_tol=5e-2, plot=False, iprint=True)
#print(tc/w_e)
print('Runtime = %g' % (time() - start))
