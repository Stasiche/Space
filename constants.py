gamma = 0.01
a = 0.01
beta = 1
D = 1
alpha = 0.001

G = 1

epsilon = 1e-10

theta = 0.3

k = 5     # количество рядов до середины
l = 5     # количество в крайнем ряду

box_size = a * k * 4
npar = int(9)

dt = 1e-6

steps_number = int(1e4)

method = 'bh+brue'
# method = 'bh'
# method = 'brute'

realization_num = 1

median = 0
