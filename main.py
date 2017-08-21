# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:40:00 2017
#
@author: Connor
"""
#import mu_zeta_chi
#import phi_matrix
import numpy as np
import tc_func as tf
from time import time
from time import strftime
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

"""
other params
"""

suffix = ''
dos_type = 'NN'
directory = './dat_files/'
date = strftime("%Y_%m_%d")

w_e = 1.0
lam_want = 0.25*tf.ttp
dom_lim = 75
num  = 15
n_lower = 0.9
n_upper = 1.0
tol = 1e-5
damp = 0.1
maxiter = 200
n = 1.0
# n_domain = np.linspace(n_lower, n_upper, num)
# tc_over_w = np.zeros(num)
# mu = np.zeros(num)
# chi_array = np.zeros([dom_lim, num])
# zeta_array = np.zeros([dom_lim, num])
# phi_array = np.zeros([dom_lim, num])
# z_array = np.zeros([dom_lim, num])

# np.savetxt((directory + dos_type + '_tc_vs_n_' + date + '_n_' +  str(n_lower) +  '_to_' +  str(n_upper) + 
#             '_steps_' +  str(num)  + '_tc_vals_lam_' + str(lam_want) + '_we_' + str(w_e) + '_dom_lim_' + 
#             str(dom_lim)  + '_tol_' +  str(tol) + '.dat'), tc_over_w)

# for c, n in enumerate(n_domain):
#    print('%g%%' % ((c+1)/num*100))
#    print('n = %g' % n_domain[c])
#    [tc_over_w[c], mu[c], chi_array[:, c], zeta_array[:, c],
#     phi_array[:, c], z_array[:,c]] = tc_solver(
#         lam_want, w_e, n, dom_lim, maxiter=maxiter,
#         mu_tol=tol, damp=damp, iprint=False)
#    tc_over_w[c] = tc_over_w[c]/w_e

#    np.savetxt((directory + dos_type + '_tc_vs_n_' + date + '_n_' +  str(n_lower) +  '_to_' +  str(n_upper) + 
#                '_steps_' +  str(num)  + '_tc_vals_lam_' + str(lam_want) + '_we_' + str(w_e) + '_dom_lim_' + 
#                str(dom_lim)  + '_tol_' +  str(tol) + '.dat'), tc_over_w)
#    np.savetxt((directory + dos_type + '_tc_vs_n_' + date + '_n_' +  str(n_lower) +  '_to_' +  str(n_upper) + 
#                '_steps_' +  str(num)  + '_n_domain_lam_' + str(lam_want) + '_we_' + str(w_e) + '_dom_lim_' + 
#                str(dom_lim)  + '_tol_' +  str(tol) + '.dat'), n_domain)
#    np.savetxt((directory + dos_type + '_tc_vs_n_' + date + '_n_' +  str(n_lower) +  '_to_' +  str(n_upper) + 
#                '_steps_' +  str(num)  + '_mu_lam_' + str(lam_want) + '_we_' + str(w_e) + '_dom_lim_' +
#                str(dom_lim)  + '_tol_' +  str(tol) + '.dat'), mu)
#    np.savetxt((directory + dos_type + '_tc_vs_n_' + date + '_n_' +  str(n_lower) +  '_to_' +  str(n_upper) + 
#                '_steps_' +  str(num)  + '_chi_lam_' + str(lam_want) + '_we_' + str(w_e) + '_dom_lim_' +
#                str(dom_lim)  + '_tol_' +  str(tol) + '.dat'), chi_array)
#    np.savetxt((directory + dos_type + '_tc_vs_n_' + date + '_n_' +  str(n_lower) +  '_to_' +  str(n_upper) + 
#                '_steps_' +  str(num)  + '_phi_lam_' + str(lam_want) + '_we_' + str(w_e) + '_dom_lim_' +
#                str(dom_lim)  + '_tol_' +  str(tol) + '.dat'), phi_array)
#    np.savetxt((directory + dos_type + '_tc_vs_n_' + date + '_n_' +  str(n_lower) +  '_to_' +  str(n_upper) + 
#                '_steps_' +  str(num)  + '_z_lam_' + str(lam_want) + '_we_' + str(w_e) + '_dom_lim_' +
#                str(dom_lim)  + '_tol_' +  str(tol) + '.dat'), z_array) 
print("Lambda %g w_e %g n %g dom_lim %g" %(lam_want, w_e, n, dom_lim))
print(tc_solver(
      lam_want, w_e, n, dom_lim, maxiter=maxiter, mu_tol = tol, damp=damp)[0]/w_e)

print('Runtime = %g' % (time() - start))
print('*************************************')
