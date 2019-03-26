from node import Node
from particle import Particle
from vector3d import Vector3d
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

from mpl_toolkits.mplot3d import Axes3D
import sys
from numpy import linalg

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


class EllipsoidTool:
    """Some stuff for playing with ellipsoids"""

    def __init__(self):
        pass

    def getMinVolEllipse(self, P=None, tolerance=0.01):
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
            M = np.diag(np.dot(QT, np.dot(linalg.inv(V), Q)))  # M the diagonal vector of an NxN matrix
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
        A = linalg.inv(
            np.dot(P.T, np.dot(np.diag(u), P)) -
            np.array([[a * b for b in center] for a in center])
        ) / d

        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0 / np.sqrt(s)

        return (center, radii, rotation)

    def getEllipsoidVolume(self, radii):
        """Calculate the volume of the blob"""
        return 4. / 3. * np.pi * radii[0] * radii[1] * radii[2]

    def plotEllipsoid(self, center, radii, rotation, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2):
        """Plot an ellipsoid"""
        make_ax = ax == None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)

        # cartesian coordinates that correspond to the spherical angles:
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rotation) + center

        if plotAxes:
            # make some purdy axes
            axes = np.array([[radii[0], 0.0, 0.0],
                             [0.0, radii[1], 0.0],
                             [0.0, 0.0, radii[2]]])
            # rotate accordingly
            for i in range(len(axes)):
                axes[i] = np.dot(axes[i], rotation)

            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                Z3 = np.linspace(-p[2], p[2], 100) + center[2]
                ax.plot(X3, Y3, Z3, color=cageColor)

        # plot ellipsoid
        ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color=cageColor, alpha=cageAlpha)

        if make_ax:
            plt.show()
            plt.close(fig)
            del fig


# make 100 random points
input_path = '../3d/clust2.xyz'
save_path = os.path.dirname(os.path.realpath(__file__)) + '/test'

if not os.path.exists(save_path):
    os.makedirs(save_path)

# P = np.reshape([random.random() * 100 for i in range(300)], (100, 3))
particles = io_xyz.read(input_path, format='Nick')
io_xyz.write(particles, save_path, 0)

# find the ellipsoid
ET = EllipsoidTool()
(center, radii, rotation) = ET.getMinVolEllipse(
    np.array([[particle.r.x, particle.r.y, particle.r.z] for particle in particles], dtype='float32'), .01)

num_sumples = 5000
k = 0

def criterier(*args):
    return True

particles = []
# for i in range(len(center)):
#     center[i] = 0
#     radii[i] = 10*(i+1)
#     for j in range(len(center)):
#         if i == j:
#             rotation[i, j] = 1
#         else:
#             rotation[i, j] = 0

while True:
    u = np.random.sample() * 2.0 * np.pi
    v = np.random.sample() * np.pi
    r1 = np.random.sample() * radii[0]
    r2 = np.random.sample() * radii[1]
    r3 = np.random.sample() * radii[2]

    tmp = [r1 * np.outer(np.cos(u), np.sin(v)), r2 * np.outer(np.sin(u), np.sin(v)), r3 * np.outer(np.ones_like(u), np.cos(v))]
    if criterier(tmp):
        [x, y, z] = np.dot([tmp[0][0][0], tmp[1][0][0], tmp[2][0][0]], rotation) + center
        particles.append(Particle(r=Vector3d(x, y, z), color=(0, 0, 0), name=len(particles)))

        k += 1
    if k % 1000 == 0:
        print(k)
    if k == num_sumples:
        break

# io_xyz.append_write(particles, save_path, 0)
io_xyz.write(particles, save_path, -999)



# particles = make_hex(constants.k, constants.l)

# for step in range(constants.steps_number):
#     print(step)
#
#     save_path = os.path.dirname(os.path.realpath(__file__)) + '/test'
#     if not os.path.exists(save_path):
#         os.makedirs(save_path)
#     with open(save_path + '/results_' + str(step) + '.xyz', 'w') as outfile:
#         outfile.write(str(len(particles)) + '\n\n')
#         for i, particle in enumerate(particles):
#             outfile.write(str(particle.r) + ' ' +
#                           str(particle.color[0]) + ' ' +
#                           str(particle.color[1]) + ' ' +
#                           str(particle.color[2]) + ' ' + '\n')
#
#     for i, particle in enumerate(particles):
#         brute_force = Vector2d(0, 0)
#         for particle2 in particles:
#             if particle != particle2:
#                 brute_force += forces.base_force(particle, particle2)
#
#         particle.force = brute_force
#
#     # print(bh.save_particles, bh.save_particles/(len(particles)*(len(particles) - 1)))
#
#     for particle in particles:
#         particle.v += particle.force * constants.dt
#         particle.r += particle.v * constants.dt
#         particle.v = Vector3d()


