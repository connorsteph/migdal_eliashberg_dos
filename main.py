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
D = 1/(emax-emin)
w_e = 16/2.5
lam_want = 2
g = np.sqrt(lam_want*w_e/2/tf.dos[np.int(tf.nee/2 + 1)])

"""
Algorithm params
"""

dom_lim = 200
p_damp = 0.3
maxiter = 50
tol = 1e-5
damp = 0.3


plt.figure()
plt.plot(np.linspace(1, 48, 10), [tc_calc(g, w_e, D, dom_lim, maxiter=maxiter,
        tol=tol, p_tol=5e-5, t_tol=5e-2, plot=False, iprint=False)
        for w_e in np.linspace(1, 48, 10)])

#tc = tc_calc(g, w_e, D,
#        dom_lim, maxiter=maxiter,
#        tol=tol, p_tol=5e-5, t_tol=5e-2, plot=False, iprint=False)
#print(tc/w_e)
print('Runtime = %g' % (time() - start))
