
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:08:35 2017

@author: Connor
"""
import tc_func as tf
from zeta_solver import zeta_solver
import numpy as np
import scipy.integrate as itg
import matplotlib.pyplot as plt

g=1.0
w_e=1.0
t0=0.1
maxiter=50
tol=1e-5
dom_lim=200
eval_lim=175
"""
begin with a guess for Tc, and then converge zeta (normalized to 1 over the interval
it was calculated on )
"""

zeta1 = zeta_solver(t0,g,w_e,dom_lim,eval_lim,maxiter,tol,iprint=False)
zeta2 = zeta_solver(t0,g,5.0,dom_lim,eval_lim,maxiter,tol,iprint=False)
zeta3 = zeta_solver(t0,3,w_e,dom_lim,eval_lim,maxiter,rel_tol,iprint=False)
zeta4 = zeta_solver(t0,5,w_e,dom_lim,eval_lim,maxiter,rel_tol,iprint=False)

plt.figure(3)
plt.hold(True)
plt.plot(tf.f_range(1,dom_lim,t0),zeta1(tf.f_range(1,dom_lim,t0)),label='zeta g=1,w_e=1')
plt.plot(tf.f_range(1,dom_lim,t0),zeta2(tf.f_range(1,dom_lim,t0)),label='zeta g=1,w_e=5')
plt.plot(tf.f_range(1,dom_lim,t0),zeta3(tf.f_range(1,dom_lim,t0)),label='zeta g=3,w_e=1')
plt.plot(tf.f_range(1,dom_lim,t0),zeta4(tf.f_range(1,dom_lim,t0)),label='zeta g=5,w_e=1')
plt.hold(False)
plt.legend(loc='best')
plt.savefig('zeta_strengths.pdf',bbox_inches='tight')
