import numpy as np
from mu_solver import mu_solver
import tc_func as tf


w_e = 0.5
t = 0.05/w_e
llam = 0.25
n = 1.0
emin = tf.e_min
emax = tf.e_max
y1 = []
y2 = []
x = []
for Nc in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
    init_zeta = [tf.freq_m(i, t) for i in range(1, Nc + 1)]
    init_chi = [tf.init_chi(w) for w in tf.freq_array(1, Nc, t)]
    dos_mu = 1/(emax - emin)
    damp = 0.3
    maxiter = 150
    tol = 1e-5
    mu, zeta, chi, z = mu_solver(
        t, llam, w_e, n, init_chi, init_zeta, Nc,
        maxiter=maxiter, damp=damp, tol=tol)

    x.append(Nc)
    y1.append(zeta[0])
    y2.append(zeta[9])
    np.savetxt('./dat_files/zeta_convergence/bcc_NNN_zeta_1_v2.dat', y1)
    np.savetxt('./dat_files/zeta_convergence/bcc_NNN_zeta_10_v2.dat', y2)
    np.savetxt('./dat_files/zeta_convergence/bcc_NNN_domain_v2.dat', x) 
    print(Nc, zeta[0], zeta[9])
