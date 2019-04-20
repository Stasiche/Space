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


def criterion(tmp, particles):
    tmp_particle = Particle(r=Vector3d(*tmp))
    for particle in particles:
        # print((tmp_particle.r - particle.r).dot(forces.base_force(tmp_particle, particle)))
        # if (tmp_particle.r - particle.r).dot(forces.base_force(tmp_particle, particle)) > -0.3:
        if abs(tmp_particle.r - particle.r) < 0.5*constants.a_0:
            return False
    return True


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


def get_and_plot_density(particles, radii, rotation, center, num=100):
    cm = Vector3d()
    particles_ = copy.deepcopy(particles)

    derotation = np.linalg.inv(rotation)

    for particle in particles_:
        particle.r.x -= center[0]
        particle.r.y -= center[1]
        particle.r.z -= center[2]

        tmp = np.dot(np.array([particle.r.x, particle.r.y, particle.r.z]), derotation)
        particle.r = Vector3d(tmp[0], tmp[1], tmp[2])

        cm += particle.r
    cm = cm * (1/len(particles_))

    for i in range(len(particles)):
        particles_[i].r.x *= 1/radii[0]
        particles_[i].r.y *= 1/radii[1]
        particles_[i].r.z *= 1/radii[2]

    lengths = sorted([abs(particle.r) for particle in particles_])
    density = []
    l = []
    # particles = sorted(particles, key=lambda particle: abs(particle.r - cm))

    for i in range(len(lengths)//num):
        density.append(num*particles_[0].mass/(4/3*np.pi*(lengths[i*num + num-1]**3 - lengths[i*num]**3)))
        # l.append(lengths[i*num] + 0.5*(lengths[i*num + num-1] - lengths[i*num]))
        l.append(lengths[i*num])

    # ig, ax1 = plt.subplots(figsize=(4, 4))
    #
    # ax1.scatter(x=l, y=density, marker='o', c='r', edgecolor='b')
    # ax1.set_title('Scatter: $x$ versus $y$')
    # ax1.set_xlabel('$x$')
    # ax1.set_ylabel('$y$')
    # plt.show()

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
            # if cur_num_particles % 1 == 0:
            #     print(cur_num_particles)
            if cur_num_particles == num//10:
                break
    return particles


input_path = '../new/3d/clusters/clust1804.xyz'
# input_path = '../3d/test/results_-999.xyz'
save_path = os.path.dirname(os.path.realpath(__file__)) + '/test_good'

if not os.path.exists(save_path):
    os.makedirs(save_path)

particles = io_xyz.read(input_path, mode='Nick')

print('Number of particles before: ', len(particles))

mean_vel = Vector3d()
kinetic_energy = 0
for particle in particles:
    mean_vel += particle.v
    kinetic_energy += particle.mass * abs(particle.v)**2 / 2
mean_vel = mean_vel * (1/len(particles))

omega = 0
for particle in particles:
    R = np.sqrt(particle.r.y ** 2 + particle.r.z ** 2)
    omega += abs(particle.v - mean_vel) / R
omega /= len(particles)

(center, radii, rotation) = getMinVolEllipse(
    np.array([[particle.r.x, particle.r.y, particle.r.z] for particle in particles], dtype='float32'))


dens, lengths = get_and_plot_density(particles, radii, rotation, center)

# for i in range(len(radii)):
    #!!!!!!!!!!!!!!
    # radii[i] = radii[i] / 5
# particles = create_sphere(radii, 2000)
# !!! del
for particle in particles:
    particle.r -= Vector3d(center[0], center[1], center[2])
io_xyz.write(particles, save_path, -1)

particles = create_sphere2(radii, lengths, mass=10)


######
(center, radii, rotation) = getMinVolEllipse(
    np.array([[particle.r.x, particle.r.y, particle.r.z] for particle in particles], dtype='float32'))

dens, lengths = get_and_plot_density(particles, radii, rotation, center, 10)


# particles = [Particle(r=Vector3d(0,0,0)),Particle(r=Vector3d(1.5*constants.a,0,0))]
# criterion([particles[0].r.x, particles[0].r.y, particles[0].r.z],[particles[1]])


print(len(particles))

for particle in particles:
    R = np.sqrt(particle.r.y ** 2 + particle.r.z ** 2)
    v_mag = omega * R
    e = Vector3d()
    e.x = 0
    e.y = particle.r.y
    e.z = particle.r.z
    tmp_v = e * (1 / abs(e)) * v_mag

    # rotate pi/2
    tmp = tmp_v.y
    tmp_v.y = tmp_v.z
    tmp_v.z = -tmp

    particle.v = tmp_v

mean_vel = Vector3d()
for particle in particles:
    mean_vel += particle.v
mean_vel *= (1/len(particles))
for particle in particles:
    particle.v -= mean_vel

print(kinetic_energy)
kinetic_energy = 0
for particle in particles:
    kinetic_energy += particle.mass * abs(particle.v)**2 / 2
print(kinetic_energy)

for step in range(constants.steps_number):
    print('s', step)
    if step % 10 == 0:
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

#
