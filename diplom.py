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

def make_univer(n_particles=int(10e6)):
    max_v = 1
    min_v = -1
    particles = [Particle(r=Vector2d(random.uniform(-box_size // 2, box_size // 2),
                                             random.uniform(-box_size // 2, box_size // 2)),
                                  v=Vector2d(random.uniform(min_v, max_v),
                                             random.uniform(min_v, max_v)))]
    # global box_size

    i = 1
    while True:
        # if i % (n_particles//100) == 0:
        #     print('box:', i // (n_particles//100) + 1)
        tmp_par = Particle(r=Vector2d(random.uniform(-box_size // 2, box_size // 2),
                                             random.uniform(-box_size // 2, box_size // 2)),
                                  v=Vector2d(random.uniform(min_v, max_v),
                                             random.uniform(min_v, max_v)))
        flag = 1
        # for particle1 in particles:
        #     tmp_e = particle1.r - tmp_par.r
        #     if abs(tmp_e) == 0:
        #         continue
        #     tmp_f = forces.base_force(particle1, tmp_par)
        #     if (tmp_f.x*tmp_e.x + tmp_f.y*tmp_e.y) < 0:
        #         flag = 0

        if flag == 1:
            particles.append(tmp_par)
            i += 1

        if i >= n_particles:
            break

        # print(i)
    return particles


def make_net(n_particles=int(10e6)):
    particles = []
    n_rows = int(math.sqrt(n_particles))
    h = constants.box_size // n_rows
    for i in range(n_particles):
        particles.append(Particle(r=Vector2d(i % n_rows * h, i // n_rows * h), v=Vector2d()))
        # print(particles[-1].r)

    return particles, h


def make_hex(k=10, l=10):
    # k количество рядов до середины
    # l количество в крайнем ряду
    particles = []
    for i in range(k):
        for j in range(i + l):
            particles.append(
                Particle(r=Vector2d(constants.a * i * math.sqrt(3) / 2,
                                    constants.a * j - constants.a * 0.5 * i),
                         v=Vector2d(), name=len(particles)))
    for i in range(k, 2 * k - 1):
        for j in range(l + k - 2 - (i - k)):
            particles.append(Particle(r=Vector2d(constants.a * i * math.sqrt(3) / 2,
                                                 constants.a * j + constants.a * 0.5 * (i - 2 * k + 2)),
                                      v=Vector2d(), name=len(particles)))
    return particles


def make_hex_2(r, center):
    particles = [Particle(r=center, v=Vector2d(), name=0)]

    for i in range(r):
        tmp_par = []
        for angle in [j * math.pi / 3 for j in range(6)]:
            tmp_par.append(Particle(r=center + Vector2d((i+1)*constants.a * math.sin(angle), (i+1)*constants.a * math.cos(angle)),
                         v=Vector2d(), name=len(particles)))

        for j, particle in enumerate(tmp_par[:-1]):
            particles.append(particle)

            # fig, ax1 = plt.subplots(figsize=(4, 4))
            # x = [particle.r.x for particle in particles]
            # y = [particle.r.y for particle in particles]
            # ax1.scatter(x=x, y=y, marker='o', c='r', edgecolor='b')
            # ax1.set_title('Scatter: $x$ versus $y$')
            # ax1.set_xlabel('$x$')
            # ax1.set_ylabel('$y$')
            # plt.show()

            e_tmp = particle.r - tmp_par[j+1].r
            e_tmp = e_tmp * (1/abs(e_tmp))
            for k in range(i):
                particles.append(Particle(r=particle.r - e_tmp * constants.a, v=Vector2d(), name=len(particles)))

                # fig, ax1 = plt.subplots(figsize=(4, 4))
                # x = [particle.r.x for particle in particles]
                # y = [particle.r.y for particle in particles]
                # ax1.scatter(x=x, y=y, marker='o', c='r', edgecolor='b')
                # ax1.set_title('Scatter: $x$ versus $y$')
                # ax1.set_xlabel('$x$')
                # ax1.set_ylabel('$y$')
                # plt.show()

        particles.append(tmp_par[-1])
        e_tmp = tmp_par[-1].r - tmp_par[0].r
        e_tmp = e_tmp * (1 / abs(e_tmp))
        for k in range(i):
            particles.append(Particle(r=particle.r - e_tmp * constants.a, v=Vector2d(), name=len(particles)))
    return particles

mean_errors = []
box_size = constants.box_size
npar = constants.npar

flag = True
list_cos_error = []
list_len_error = []

particles = make_hex(constants.k, constants.l)
# particles = make_hex_2(3, Vector2d())

fig, ax1 = plt.subplots(figsize=(4, 4))
x = [particle.r.x for particle in particles]
y = [particle.r.y for particle in particles]
ax1.scatter(x=x, y=y, marker='o', c='r', edgecolor='b')
ax1.set_title('Scatter: $x$ versus $y$')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$y$')
plt.show()

constants.median = len(particles) // 2

# eps = 0.5
# for i in range(-100, 11):
#     if constants.a + i * eps != 0:
#         f = forces.base_force(Particle(r=Vector2d(), name='l'),
#                               Particle(r=Vector2d(constants.a + i * eps, 0), name='r'))
#         print(np.round((constants.a + i * eps)/constants.a, 3), f)

for step in range(constants.steps_number):
    print(step)

    save_path = os.path.dirname(os.path.realpath(__file__)) + '/test'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    with open(save_path + '/results_' + str(step) + '.xyz', 'w') as outfile:
        outfile.write(str(len(particles)) + '\n\n')
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.color[0]) + ' ' +
                          str(particle.color[1]) + ' ' +
                          str(particle.color[2]) + ' ' + '\n')

    # fig, ax1 = plt.subplots(figsize=(4, 4))
    # x = [particle.r.x for particle in particles]
    # y = [particle.r.y for particle in particles]
    # ax1.scatter(x=x, y=y, marker='o', c='r', edgecolor='b')
    # ax1.set_title('Scatter: $x$ versus $y$')
    # ax1.set_xlabel('$x$')
    # ax1.set_ylabel('$y$')
    # plt.show()

    # mean_error = 0
    # n1 = Node([-box_size / 2, -box_size / 2, box_size / 2, box_size / 2], node_type=2)
    # bh = BHtree(n1)
    #
    # for particle in particles:
    #     bh.insert_particle(particle, bh.root)
    #
    # if step == 0:
    #     plotting.plot_main(particles, bh)

    for i, particle in enumerate(particles):
        # if i == constants.median:
        #     q = 1
        # bh_force = bh.calculate_force(particle, bh.root)

        brute_force = Vector2d(0, 0)
        for particle2 in particles:
            if particle != particle2:
                brute_force += forces.base_force(particle, particle2)

        # bh_force_abs = abs(bh_force)
        # brute_force_abs = abs(brute_force)
        #
        # len_error = (bh_force_abs - brute_force_abs) / (brute_force_abs + constants.epsilon)
        #
        # # print(bh_force, '|', brute_force)
        # tmp1 = bh_force.x/(brute_force.x + constants.epsilon)
        # tmp2 = bh_force.y/(brute_force.y + constants.epsilon)
        # eps = 1e-5
        # if ((abs(tmp1-1) > 0.1) or (abs(tmp2-1) > 0.1)) \
        #         and not((abs(bh_force.x) < eps) or (abs(bh_force.y) < eps) or (abs(brute_force.x) < eps) or (abs(brute_force.x) < eps)):
        #     print(i)
        #     print(bh_force, '|', brute_force)
        #     print(bh_force.x/(brute_force.x + constants.epsilon), '|', bh_force.y/(brute_force.y + constants.epsilon))
        #     print('____')

        # particle.force = bh_force
        particle.force = brute_force

    # print(bh.save_particles, bh.save_particles/(len(particles)*(len(particles) - 1)))

    for particle in particles:
        particle.v += particle.force * constants.dt
        particle.r += particle.v * constants.dt
        particle.v = Vector2d()

