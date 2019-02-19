import matplotlib.pyplot as plt
# import seaborn as sns
from vector2d import Vector2d
from particle import Particle
from node import Node
import random


def make_univer(n_particles=int(10e6)):
    particles = []
    # global box_size
    max_v = 1
    min_v = -1
    for i in range(n_particles):
        if i % (n_particles//100) == 0:
            print(i // (n_particles//100) + 1)
        particles.append(Particle(r=Vector2d(random.randint(-box_size // 2, box_size // 2),
                                             random.randint(-box_size // 2, box_size // 2)),
                                  v=Vector2d(random.randint(min_v, max_v),
                                             random.randint(min_v, max_v))))

    return particles


box_size = 10
npar = int(10e6)
particles = make_univer(npar)
mvx = 0
mvy = 0
mrx = 0
mry = 0
for el in particles:
    mvx += el.v.x
    mvy += el.v.y
    mrx += el.r.x
    mry += el.r.y

mvx /= npar
mvy /= npar
mrx /= npar
mry /= npar

print((mvx**2 + mvy**2)**0.5, (mrx**2 + mry**2)**0.5)

