# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:40:00 2017

@author: Connor
"""
#import mu_zeta_chi
import numpy as np
import tc_func as tf
from time import time
import matplotlib.pyplot as plt
from mu_solver import mu_solver
from tc_solver import tc_solver
start = time()
emin = tf.e_min
emax = tf.e_max
dos_avg = tf.dos_avg
"""
Problem params
"""
ttp = tf.ttp
D = (emax-emin)
w_e = 1*tf.ttp
lam_want = 2*tf.ttp
n = 0.15
"""
Algorithm params
"""
dom_lim = 150
maxiter = 150
tol = 1e-3
damp = 0.2
#0.17 for const_dos
#0.2 for bcc_dos
plt.close('all')
#********************************************************
#Nc = dom_lim+30
#damp = 0.1
#mu = 5.
#dos = tf.dos
#t = 1.0*w_e
#init_chi = [tf.init_chi(w) for w in tf.freq_array(1, dom_lim, t)]
#init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, dom_lim, t)]
#zeta, chi, iterations = mu_zeta_chi.simult_mu_func(t,lam_want,w_e,mu,tf.dos_avg,emin,emax,tf.dee,damp,
#                                           maxiter,tol,tf.dos,init_zeta,init_chi)
#plt.figure()
#plt.plot(zeta)
#plt.figure()
#plt.plot(chi)
#print(iterations)
#************************************************************
num = 30
n_domain = np.linspace(0.15, 1.0, num)
y = np.empty(num)
for c, i in enumerate(n_domain):
    print('%g%%' % ((c+1)/num*100))
    y[c] = tc_solver(
        lam_want, w_e, i, dom_lim, maxiter=maxiter,
        tol=tol, p_tol=tol, iprint=False)[0]
    y[c] = y[c]/w_e
    np.savetxt('./dat_files/bcc_NN_t_vs_n_2017_06_30_n_0_15_to_1_steps_30.dat', y)
np.savetxt('./dat_files/bcc_NN_t_vs_n_2017_06_30_n_0_15_to_1_steps_30.dat', y)
plt.figure()
plt.grid(True)
plt.plot(n_domain, y)
plt.xlabel(r'$n$', fontsize=14)
plt.xlim([0,1.0])
plt.ylim([0,0.5])
plt.ylabel(r'$\frac{T_c}{\omega_E}$', fontsize=14)
plt.title('Tc versus filling for BCC DOS (NN)\n w_e = %g, lam = %g' % (w_e, lam_want))
plt.savefig('./plots/tc_vs_n_bcc_dos_ttx_0.pdf', dpi=150, bbox_inches='tight')
#************************************************************
#tc, mu, chi, zeta, phi, z = tc_solver(
#        lam_want, w_e, n, dom_lim, maxiter=maxiter,
#        tol=tol, p_tol=tol, iprint=False)
#
#plt.figure()
#plt.plot(phi)
#plt.title('phi')
#plt.figure()
#plt.plot(chi)
#plt.title('chi')
#plt.figure()
#plt.plot(zeta)
#plt.title('zeta')
#print('\n********************************\nLambda = %g' %lam_want)
#print('n is %g' % n)
#print('w_e is %g' % w_e)
#print('Mu is %g' % mu)
#print('Tc/we is %g' % (tc/w_e))
#******************************************
#tc = 0.9*w_e
#y=[]
#for damp in np.arange(0,0.99,0.05):
#    init_chi = [tf.init_chi(w) for w in tf.freq_array(1, dom_lim, tc)]
#    init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, dom_lim, tc)]
#    start2 = time()
#    mu, zeta, chi, z = mu_solver(tc, lam_want, w_e, n, init_chi, init_zeta,
#                          dom_lim, maxiter=maxiter, tol=tol, damp=damp)
#    y.append(time()-start2)
#plt.figure()
#plt.grid(True)
#plt.plot(np.arange(0,0.99,0.05),y)

#init_chi = [tf.init_chi(w) for w in tf.freq_array(1, dom_lim, tc)]
#init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, dom_lim, tc)]
#mu, zeta, chi, z = mu_solver(tc, lam_want, w_e, n, init_chi, init_zeta,
#                          dom_lim, maxiter=maxiter, tol=tol, damp=damp)
#print(mu)
#plt.figure()
#plt.plot(z)
#plt.figure()
#plt.plot(chi)
#plt.ylim([-3,3])
print('Runtime = %g' % (time() - start))
