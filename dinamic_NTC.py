from node import Node
from particle import Particle
from vector2d import Vector2d
from bh_tree import BHtree
from profiler import Profiler
import random
import math
import forces
import constants

import networkx as nx
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


mean_errors = []
box_size = constants.box_size
npar = constants.npar

flag = True
list_cos_error = []
list_len_error = []
for num_realization in range(constants.realization_num):
    if flag:
        particles = make_hex(10, 10)
        # particles = make_univer(npar)
        # particles, h = make_net(npar)
    else:
        # Шестиугольник
        # particles = []
        # for i in range(3):
        #     for j in range(i+3):
        #         particles.append(
        #             Particle(r=Vector2d(constants.a * i * math.sqrt(3)/2,
        #                                 constants.a * j - constants.a * 0.5 * i),
        #                      v=Vector2d(), name=len(particles)))
        # for i in range(3, 5):
        #     for j in range(7-i):
        #         particles.append(Particle(r=Vector2d(constants.a * i * math.sqrt(3)/2,
        #                                              constants.a * j + constants.a * 0.5 * (i-4)),
        #                                   v=Vector2d(), name=len(particles)))
        ###################################

        particles = []
        k = 10              # количество рядов до середины
        l = 10               # количество в крайнем ряду
        for i in range(k):
            for j in range(i + l):
                particles.append(
                    Particle(r=Vector2d(constants.a * i * math.sqrt(3) / 2,
                                        constants.a * j - constants.a * 0.5 * i),
                             v=Vector2d(), name=len(particles)))
        for i in range(k, 2*k - 1):
            for j in range(l+k - 2 - (i - k)):
                particles.append(Particle(r=Vector2d(constants.a * i * math.sqrt(3) / 2,
                                                     constants.a * j + constants.a * 0.5 * (i - 2*k + 2)),
                                          v=Vector2d(), name=len(particles)))
    # fig, (ax1, ax2) = plt.subplots(
    #     nrows=1, ncols=2,
    #     figsize=(8, 4)
    # )
    #
    fig, ax1 = plt.subplots(figsize=(4, 4))
    x = [particle.r.x for particle in particles]
    y = [particle.r.y for particle in particles]
    ax1.scatter(x=x, y=y, marker='o', c='r', edgecolor='b')
    ax1.set_title('Scatter: $x$ versus $y$')
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$')

    plt.show()

    mean_error = 0
    n1 = Node([-box_size / 2, -box_size / 2, box_size / 2, box_size / 2], node_type=2)
    bh = BHtree(n1)

    for particle in particles:
        bh.insert_particle(particle, bh.root)

    for i, particle in enumerate(particles):
        bh_force = bh.calculate_force(particle, bh.root)
        # particle.force = bh_force

        brute_force = Vector2d(0, 0)
        for particle2 in particles:
            if particle != particle2:
                brute_force += forces.base_force(particle, particle2)
                # brute_force += forces.gravity(particle, particle2)

        particle.force = bh_force

        # print(bh_force, '\n', brute_force, '\n', abs(bh_force - brute_force), '\n', '______')
        # print(abs(bh_force - brute_force))
        bh_force_abs = abs(bh_force)
        brute_force_abs = abs(brute_force)

        cos_error = bh_force.dot(brute_force) / (bh_force_abs * brute_force_abs + constants.epsilon)
        len_error = (bh_force_abs - brute_force_abs) / (brute_force_abs + constants.epsilon)

        list_cos_error.append(cos_error)
        list_len_error.append(len_error)

        print(i)
        # print(cos_error, len_error)
        # print(len_error, '|', bh_force, '|', brute_force)
        print(bh_force.x/brute_force.x, '|', bh_force.y/brute_force.y)
        print('____')
        # print(abs(bh_force - brute_force))
        mean_error += abs(bh_force - brute_force)
    mean_error /= len(particles)
    mean_errors.append(mean_error)
    # if mean_error > 1:
    #     flag = 0

    # print(num_realization, 'mean error: ' + str(mean_error))
    # print(num_realization, 'mean error: ' + str(mean_error))

# print('cos', min(list_cos_error))
print('len', min(list_len_error), max(list_len_error))
# print('h', h)

