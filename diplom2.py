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
import copy
import os

# import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def make_rotation_matrix(angle):
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


def getMinVolEllipse(P=None, tolerance=0.01):
    """ Find the minimum volume ellipsoid which holds all the points

    Based on work by Nima Moshtagh
    http://www.mathworks.com/matlabcentral/fileexchange/9542
    and also by looking at:
    http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
    Which is based on the first reference anyway!

    Here, P is a numpy array of N dimensional points like this:
    P = [[x,y,z,...], <-- one point per line
         [x,y,z,...],
         [x,y,z,...]]

    Returns:
    (center, radii, rotation)

    """
    (N, d) = np.shape(P)
    d = float(d)

    # Q will be our working array
    Q = np.vstack([np.copy(P.T), np.ones(N)])
    QT = Q.T

    # initializations
    err = 1.0 + tolerance
    u = (1.0 / N) * np.ones(N)

    # Khachiyan Algorithm
    while err > tolerance:
        V = np.dot(Q, np.dot(np.diag(u), QT))
        M = np.diag(np.dot(QT, np.dot(np.linalg.inv(V), Q)))  # M the diagonal vector of an NxN matrix
        j = np.argmax(M)
        maximum = M[j]
        step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
        new_u = (1.0 - step_size) * u
        new_u[j] += step_size
        err = np.linalg.norm(new_u - u)
        u = new_u

    # center of the ellipse
    center = np.dot(P.T, u)

    # the A matrix for the ellipse
    A = np.linalg.inv(
        np.dot(P.T, np.dot(np.diag(u), P)) -
        np.array([[a * b for b in center] for a in center])
    ) / d

    # Get the values we'd like to return
    U, s, rotation = np.linalg.svd(A)
    radii = 1.0 / np.sqrt(s)

    return (center, radii, rotation)


def create_sphere(radii, num_sumples=2000):
    cur_num_particles = 0
    attempts = 100
    particles = []
    while True:
        u = np.random.uniform(0, 2*np.pi)
        # v = np.random.uniform(0, np.pi)

        # cosu = np.random.uniform(-1, 1)
        cosv = np.random.uniform(-1, 1)
        tmp_z = np.random.uniform(0, 1)

        v = np.arccos(cosv)

        mul = tmp_z**(1/3)
        r1 = radii[0]*mul
        r2 = radii[1]*mul
        r3 = radii[2]*mul

        # mul = tmp_z**(1/2)
        # r1 = radii[0]*mul
        # r2 = radii[1]*mul
        # r3 = radii[2]*mul


        # lam = 0.1
        # k = 0.00001
        # r1 = radii[0]*(1-np.exp(tmp_z/lam)**k)
        # r2 = radii[1]*(1-np.exp(tmp_z/lam)**k)
        # r3 = radii[2]*(1-np.exp(tmp_z/lam)**k)

        # k = 1
        # mul = -1/np.log(k - k*0.99999) * np.log(k - k*tmp_z)
        # r1 = radii[0]*mul
        # r2 = radii[1]*mul
        # r3 = radii[2]*mul

        # k = 1
        # mul = (-1/(1+np.exp(-10*(tmp_z-0.5)))+1)**(1/3)
        # r1 = radii[0]*mul
        # r2 = radii[1]*mul
        # r3 = radii[2]*mul

        tmp = [r1 * np.cos(u) * np.sin(v), r2 * np.sin(u) * np.sin(v), r3 * np.cos(v)]
        if criterion(tmp, particles):
            # [x, y, z] = np.dot([tmp[0][0][0], tmp[1][0][0], tmp[2][0][0]], rotation) + center
            [x, y, z] = [tmp[0], tmp[1], tmp[2]]
            particles.append(Particle(r=Vector3d(x, y, z), color=(0, 0, 0), name=len(particles)))

            cur_num_particles += 1
        if cur_num_particles % 1 == 0:
            print(cur_num_particles)
        if cur_num_particles == num_sumples:
            break
    return particles


def get_and_plot_density(particles, radii, rotation, center, num=5):
    cm = Vector2d()
    particles_ = copy.deepcopy(particles)

    derotation = np.linalg.inv(rotation)

    for particle in particles_:
        particle.r.x -= center[0]
        particle.r.y -= center[1]

        tmp = np.dot(np.array([particle.r.x, particle.r.y]), derotation)
        particle.r = Vector2d(tmp[0], tmp[1])

        cm += particle.r
    cm = cm * (1/len(particles_))

    for i in range(len(particles)):
        particles_[i].r.x *= 1/radii[0]
        particles_[i].r.y *= 1/radii[1]

    lengths = sorted([abs(particle.r) for particle in particles_])
    density = []
    l = []
    # particles = sorted(particles, key=lambda particle: abs(particle.r - cm))

    for i in range(len(lengths)//num):
        density.append(num*particles_[0].mass/(np.pi*(lengths[i*num + num-1]**2 - lengths[i*num]**2)))
        # l.append(lengths[i*num] + 0.5*(lengths[i*num + num-1] - lengths[i*num]))
        l.append(lengths[i*num])

    ig, ax1 = plt.subplots(figsize=(4, 4))

    ax1.scatter(x=l, y=density, marker='o', c='r', edgecolor='b')
    ax1.set_title('Scatter: $x$ versus $y$')
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$')
    plt.show()

    return density, lengths


def create_sphere2(radii, lengths, num_sumples=2000, mass=1):
    attempts = 100
    # particles = [Particle(r=Vector3d(), mass=10)]
    particles = []
    num = 100
    for i in range(0, len(lengths) // num):
        cur_num_particles = 0
        while True:
            u = np.random.uniform(0, 2*np.pi)
            cosv = np.random.uniform(-1, 1)
            tmp_z = np.random.uniform(0, 1)
            v = np.arccos(cosv)

            mul = tmp_z**(1/3)
            # mul = (mul * lengths[i*num + num-1]) + (lengths[i*num + num-1] - lengths[i*num])
            mul = lengths[i*num] + mul*(lengths[i*num + num-1] - lengths[i*num])
            r1 = radii[0]*mul
            r2 = radii[1]*mul
            r3 = radii[2]*mul

            tmp = [r1 * np.cos(u) * np.sin(v), r2 * np.sin(u) * np.sin(v), r3 * np.cos(v)]
            if criterion(tmp, particles):
                [x, y, z] = [tmp[0], tmp[1], tmp[2]]
                particles.append(Particle(r=Vector3d(x, y, z), mass=mass))

                cur_num_particles += 1
            if cur_num_particles % 1 == 0:
                print(cur_num_particles)
            if cur_num_particles == num//10:
                break
    return particles


# particles = io_xyz.read('../Space/test/results_22350.xyz')
# particles = io_xyz.read('../Space/test/results_1.xyz')
particles = io_xyz.read('../Space/test/results_4000.xyz')

(center, radii, rotation) = getMinVolEllipse(
    np.array([[particle.r.x, particle.r.y] for particle in particles], dtype='float32'))

dens, lengths = get_and_plot_density(particles, radii, rotation, center)

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
# Nice!
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


