from particle import Particle
from vector3d import Vector3d
import os
from struct import *
import numpy as np


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


def save_as_init(particles):
    with open('!init_configuration.txt', 'w') as outfile:
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.v) + ' ' +
                          str(particle.mass) + ' ' +
                          str(particle.a) + '\n\r')


path_system = os.path.dirname(os.path.realpath(__file__)) + "/out.xyz"
system = read(path_system, "out")
for particle in system:
    particle.a = 1
save_as_init(system)
