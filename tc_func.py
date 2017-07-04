# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 14:59:33 2017

@author: Connor
"""
import numpy as np

"""
important to have a data file with ttx and ttp data in the two column format
 (use the compiled program with BCC config file)
"""

data = np.loadtxt('./dat_files/ttx_0d0_nx_500_dee_5d-2.dat', skiprows=1, usecols=[0, 1])
dos_vals = data[1:, :]
[ttp, ttx] = data[0, :]
dee = dos_vals[1, 0] - dos_vals[0, 0]
nee = np.size(dos_vals[:, 0])
e_min = -8/ttp-6*abs(ttx)/ttp
e_max = 8/ttp+6*abs(ttx)/ttp
dos_avg = 1/(e_max-e_min)
cusp = 6*ttx-4*ttx**3/ttp**2
dos_domain = dos_vals[:, 0]
dos = dos_vals[:, 1]
#dos = [1/(e_max-e_min) for i in range(nee)]


def init_phi(w):
    return 1/w


def init_chi(w):
    return 1/(5+w**2)


def init_zeta(w):
    return w


def interpolater(f_domain, f_range):
    """
    given a set of x_i, f(x_i), returns an interpolated function
    f(x) defined over x in {x_i_min -> x_i_max}
    """
    def func(x):
        if isinstance(x, (list, tuple, np.ndarray)):
            vals = []
            for i in x:
                l_index = binsearch(f_domain, i)
                if l_index+1 == np.size(f_domain):
                    return
                vals.append(f_range[l_index]
                            + (i-f_domain[l_index])*(f_range[l_index+1]
                            - f_domain[l_index])/(f_range[l_index+1]
                            - f_range[l_index]))
            return vals
        else:
            l_index = binsearch(f_domain,x)
            return (f_range[l_index]
                    +(x-f_domain[l_index])*(f_range[l_index+1]
                    -f_range[l_index])/(f_domain[l_index+1]-f_domain[l_index]))
    return func


def binsearch(vec, val):
    """
    returns lower bracketing index for the closest value in vec to val
    """
    m = np.size(vec)
    ia = 0
    ib = m-1

    while (ib-ia) > 1:
        im = (ia + ib)//2
        if vec[im] < val:
            ia = im
        else:
            ib = im
    return ia


#def dos(e):
#    # interpolates the DOS from a given file
##    if isinstance(e, (list, tuple, np.ndarray)):
##        vals = []
##        for i in e:
##            l_index = binsearch(dos_vals[:, 0], i)
##            if l_index+1 == dos_vals[:, 0].size:
##                return
##            vals.append(dos_vals[l_index, 1]
##                        + (i-dos_vals[l_index, 0])*(dos_vals[l_index+1, 1]
##                        - dos_vals[l_index, 1])/(dos_vals[l_index+1, 0]
##                        - dos_vals[l_index, 0]))
##        return vals
##    else:
##        l_index=binsearch(dos_vals[:,0],e)
##        return (dos_vals[l_index,1]
##                +(e-dos_vals[l_index,0])*(dos_vals[l_index+1,1]
##                -dos_vals[l_index, 1])/(dos_vals[l_index + 1, 0]
##                - dos_vals[l_index, 0]))
#    return 1/(e_max-e_min)


def freq_m(n, t):
    """
    returns the fermion matsubara frequency corresponding to integer n
    """
    return np.pi*t*(2*n-1)


def matsu_index(w, t):
    """
    returns the integer n which gives the fermion Matsubara freq.
    w_n=pi*t(2*n-1)
    """
    return round(1/2*(w/np.pi/t+1))


def freq_array(lower, upper, t):
    """
    returns matsubara frequencies m= lower to upper
    """
    return [freq_m(n, t) for n in m_array(lower, upper)]


def m_array(lower, upper):
    """
    returns corresponding m for f_range
    """

    return range(lower, upper+1, 1)


def f_compare(v1, v2):
    """
    returns the average absolute value of the difference between
    two vector's components
    """
    diff = 0.0
    assert(np.size(v1) == np.size(v2))
    for i in range(np.size(v1)):
        diff += abs((v1[i]/v2[i])-1)
#        diff += abs(v1[i]-v2[i])
    return diff/np.size(v1)