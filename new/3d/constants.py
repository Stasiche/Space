import numpy as np

gamma_0 = 1
a_0 = 1
m_0 = 1
beta_crit = 2 * np.sqrt(11) * np.sqrt(a_0 ** 3 / (gamma_0 * m_0))

beta_ratio = 0.1

beta = beta_crit * beta_ratio
D = 1
alpha = 0.001
a_cut_0 = 2.5 * a_0
a_cut_sqr_0 = a_cut_0 * a_cut_0

G = 1

epsilon = 1e-10

theta = 0.3

k = 5     # количество рядов до середины
l = 5     # количество в крайнем ряду

box_size = a_0 * k * 4
npar = int(9)

dt = 1e-6

steps_number = int(1e6)

method = 'bh+brue'
# method = 'bh'
# method = 'brute'

realization_num = 1

median = 0
