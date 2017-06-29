# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:40:00 2017

@author: Connor
"""
import mu_zeta_chi
import numpy as np
import tc_func as tf
from time import time
import matplotlib.pyplot as plt
from mu_solver import mu_solver
from tc_solver import tc_solver
start = time()
emin = tf.e_min
emax = tf.e_max
"""
Problem params
"""
ttp = tf.ttp
D = (emax-emin)
w_e = D/2.5*tf.ttp
lam_want = 2*tf.ttp
n = 1.0
print('Lambda = %g' %lam_want)
"""
Algorithm params
"""

dom_lim = 60
maxiter = 150
tol = 1e-3

plt.close('all')
#********************************************************
#Nc = dom_lim+30
#damp = 0.1
#mu = 5.
#dos = tf.dos
#dos_mu = tf.interpolater(tf.dos_domain, tf.dos)(mu)
#g = np.sqrt(lam_want*w_e/2/dos_mu)
#t = 1.0*w_e
#init_chi = [tf.init_chi(w) for w in tf.freq_array(1, dom_lim, t)]
#init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, dom_lim, t)]
#zeta, chi, iterations = mu_zeta_chi.simult_mu_func(t,g,w_e,mu,dos_mu,emin,emax,tf.dee,damp,
#                                           maxiter,tol,tf.dos,init_zeta,init_chi)
#plt.figure()
#plt.plot(zeta)
#plt.figure()
#plt.plot(chi)
#print(iterations)
#************************************************************
tc, mu, chi, zeta, phi = tc_solver(
        lam_want, w_e, n, dom_lim, maxiter=maxiter,
        tol=tol, p_tol=tol, iprint=False)

plt.figure()
plt.plot(phi)
plt.title('phi')
plt.figure()
plt.plot(chi)
plt.title('chi')
plt.figure()
plt.plot(zeta)
plt.title('zeta')
print('n is %g' % n)
print('w_e is %g' % w_e)
print('Mu is %g' % mu)
print('Tc/we is %g' % (tc/w_e))

#******************************************
#tc = 0.01*w_e
#init_chi = [tf.init_chi(w) for w in tf.freq_array(1, dom_lim, tc)]
#init_zeta = [tf.init_zeta(w) for w in tf.freq_array(1, dom_lim, tc)]
#mu, zeta, chi = mu_solver(tc, lam_want, w_e, n, init_chi, init_zeta,
#                          dom_lim, maxiter=maxiter, tol=tol)
#print(mu)
#plt.figure()
#plt.plot(zeta)
#plt.figure()
#plt.plot(chi)
print('Runtime = %g' % (time() - start))
