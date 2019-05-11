from particle import Particle
from vector3d import Vector3d
import os
from struct import *


def read_a3r(read_path):
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


path_system = os.path.dirname(os.path.realpath(__file__)) + "/r000_s002_000000000_config.a3r"
system = read_a3r(path_system)
for particle in system:
    particle.a = 1

save_as_init(system)
