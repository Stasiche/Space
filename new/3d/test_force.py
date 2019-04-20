import matplotlib.pyplot as plt
from vector3d import Vector3d
from particle import Particle
import numpy as np
import io_xyz
import forces
import constants
import random
import copy
import os


save_path = os.path.dirname(os.path.realpath(__file__)) + '/test_force'
if not os.path.exists(save_path):
    os.makedirs(save_path)

particles = [Particle(r=Vector3d(0, 0, 0)), Particle(r=Vector3d(constants.a_0*1.00000000000001, 0, 0))]

for step in range(constants.steps_number):
    print('s', step)
    if step % 1000 == 0:
        io_xyz.write(particles, save_path, step)

    for particle in particles:
        particle.force = Vector3d()

    for i, p1 in enumerate(particles):
        for p2 in particles[i + 1:]:
            tmp_force = forces.base_force(p1, p2)
            p1.force += tmp_force
            p2.force += -tmp_force

    for particle in particles:
        particle.v += particle.force * constants.dt
        particle.r += particle.v * constants.dt
