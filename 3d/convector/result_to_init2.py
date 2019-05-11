from particle import Particle
from vector3d import Vector3d
import os
from struct import *
import numpy as np


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


def read(read_path, mode="classic"):
    particles = []
    with open(read_path, 'r') as inputfile:
        s_tmp = inputfile.read()
        for el in s_tmp.split('\n')[2:]:
            if el != '':
                if mode == "classic":
                    rx, ry, rz, vx, vy, vz = el.split()
                    particles.append(Particle(r=Vector3d(float(rx), float(ry), float(rz)),
                                              v=Vector3d(float(vx), float(vy), float(vz)), mass=2.5))
                elif mode == "clust":
                    с, rx, ry, rz, vx, vy, vz, c = el.split()
                    particles.append(Particle(r=Vector3d(float(rx), float(ry), float(rz)),
                                              v=Vector3d(float(vx), float(vy), float(vz)), mass=1))
                elif mode == "out":
                    с, rx, ry, rz, vx, vy, vz, c, c, c = el.split()
                    particles.append(Particle(r=Vector3d(float(rx), float(ry), float(rz)),
                                              v=Vector3d(float(vx), float(vy), float(vz)), mass=1))
    return particles


def read_a3r(read_path, mode=None):
    particles = []
    with open(read_path, 'rb') as inputfile:
        byte = inputfile.read(4)
        print(unpack('4s', byte))

        byte = inputfile.read(4)
        num_particles = unpack('i', byte)[0]
        print(num_particles)

        byte = inputfile.read(4)
        data_start = unpack('i', byte)[0]
        print(data_start)

        byte = inputfile.read(10)
        print(unpack('10s', byte))

        byte = inputfile.read(8)
        print(unpack('d', byte))

        byte = inputfile.read(4)
        print(unpack('i', byte))

        byte = inputfile.read(6)
        print(unpack('6s', byte))

        print("data: ")
        for i in range(num_particles):
            byte = inputfile.read(4)
            rx = unpack('f', byte)[0]

            byte = inputfile.read(4)
            ry = unpack('f', byte)[0]

            byte = inputfile.read(4)
            rz = unpack('f', byte)[0]

            particles.append(Particle(r=Vector3d(rx, ry, rz), name=len(particles)))
            print(i, rx, ry, rz)

    return particles


def save_as_init(particles):
    with open('!init_configuration.txt', 'w') as outfile:
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.v) + ' ' +
                          str(particle.mass) + ' ' +
                          str(particle.a) + '\n\r')


def write(particles, name):
    with open(name + '.xyz', 'w') as outfile:
        outfile.write(str(len(particles)) + '\n\n')
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.v) + '\n')


def read_info(path):
    with open(path + '/info.txt', 'r') as infile:
        a = float(infile.readline().split()[1])
        m = float(infile.readline().split()[1])
    return a, m


path_clusters = os.path.dirname(os.path.realpath(__file__)) + "/clusters/"
path_system = os.path.dirname(os.path.realpath(__file__)) + "/out.xyz"

for i in os.walk(path_clusters):
    if len(i[1]) == 0:
        a, m = read_info(i[0])
        print(a, m)