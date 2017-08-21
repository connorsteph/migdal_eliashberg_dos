# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 16:47:31 2017

@author: Connor
"""
import n_sum
from dot_product import vanilla_dot_product
from omp_dot_product import omp_dot_product
from time import time
from numpy.random import rand

a = rand(100000000)
b = rand(100000000)
start = time()
c = vanilla_dot_product(a, b)
print(time() - start)
start = time()
d = vanilla_dot_product(a, b)
print(time() - start)