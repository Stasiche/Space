from node import Node
from particle import Particle
from vector2d import Vector2d
from bh_tree import BHtree
from profiler import Profiler
import random
import io_xyz
import math
import forces
import constants
import plotting
import os

# import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def make_rotation_matrix(angle):
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

# particles = io_xyz.read('../Space/test/results_22350.xyz')
# particles = io_xyz.read('../Space/test/results_1.xyz')
particles = io_xyz.read('../Space/test/results_4000.xyz')

rot = make_rotation_matrix(np.pi/2)

cm = Vector2d()
for particle in particles:
    cm += particle.r
cm = cm * (1/len(particles))

# delta = 100
# for particle in particles:
#     # particle.v = Vector2d(10 * np.random.sample(), 10 * np.random.sample())
#     r = particle.r - cm
#     r = np.dot(rot, [r.x, r.y])
#     x = np.random.sample() * 2 * delta + (r[0] - delta)
#     y = np.random.sample() * 2 * delta + (r[1] - delta)
#     particle.v = Vector2d(x, y) * 1
#####
# for i, particle in enumerate(particles):
#     brute_force = Vector2d(0, 0)
#     for particle2 in particles:
#         if particle != particle2:
#             brute_force += forces.base_force(particle, particle2)
#
#     r = particle.r - cm
#     r = np.dot(rot, [r.x, r.y])
#     r = Vector2d(r[0], r[1])
#     tmp_v = r * (1/abs(r)) * abs((r.dot(brute_force)/abs(r)))
#     particle.v = tmp_v * 0.1
########
for i, particle in enumerate(particles):
    brute_force = Vector2d(0, 0)
    for particle2 in particles:
        if particle != particle2:
            brute_force += forces.base_force(particle, particle2)

    r = particle.r - cm
    r = np.dot(rot, [r.x, r.y])
    r = Vector2d(r[0], r[1])
    # tmp_v = r * (1/abs(r)) * np.sqrt(abs((r.dot(brute_force)))/particle.mass)
    tmp_v = r * (1/abs(r)) * np.sqrt(abs(brute_force) * abs(r)/particle.mass)
    angle = np.pi/12
    # particle.v = tmp_v.rotate(np.random.sample() * 2 * angle - angle)
    particle.v = tmp_v.rotate(np.random.sample() * angle)
########
# for i, particle in enumerate(particles):
#     brute_force = Vector2d(0, 0)
#     for particle2 in particles:
#         if particle != particle2:
#             brute_force += forces.base_force(particle, particle2)
#
#     r = particle.r - cm
#     tmp_v = r * (1/abs(r)) * abs((r.dot(brute_force)/abs(r)))
#     particle.v = tmp_v * 0.1
######


tmp_v = Vector2d()
for particle in particles:
    tmp_v += particle.v
tmp_v = tmp_v * (1/len(particles))

for particle in particles:
    particle.v -= tmp_v


for step in range(constants.steps_number):
    print(step)

    save_path = os.path.dirname(os.path.realpath(__file__)) + '/test_velocity'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    with open(save_path + '/results_' + str(step) + '.xyz', 'w') as outfile:
        outfile.write(str(len(particles)) + '\n\n')
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.v) + '\n')

    for i, particle in enumerate(particles):
        brute_force = Vector2d(0, 0)
        for particle2 in particles:
            if particle != particle2:
                brute_force += forces.base_force(particle, particle2)

        particle.force = brute_force

    for particle in particles:
        particle.v += particle.force * constants.dt
        particle.r += particle.v * constants.dt


