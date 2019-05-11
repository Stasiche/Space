import matplotlib.pyplot as plt
from vector3d import Vector3d
from particle import Particle
import numpy as np
import io_xyz
import forces_a
import constants
import random
import copy
import os


def get_omega(particle, rot_axis):
    tmp_pro = np.dot([particle.r.x, particle.r.y, particle.r.z], [rot_axis.x, rot_axis.y, rot_axis.z])
    R = abs(particle.r - rot_axis * tmp_pro)
    return abs(particle.v) / R


# TODO rename
def get_vel(particle, rot_axis, omega):
    tmp_pro = np.dot([particle.r.x, particle.r.y, particle.r.z], [rot_axis.x, rot_axis.y, rot_axis.z])
    R = abs(particle.r - rot_axis * tmp_pro)
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

    return tmp_v


def criterion(tmp, particles):
    # return True
    tmp_particle = Particle(r=Vector3d(*tmp))
    for particle in particles:
        # print((tmp_particle.r - particle.r).dot(forces.base_force(tmp_particle, particle)))
        # if (tmp_particle.r - particle.r).dot(forces.base_force(tmp_particle, particle)) > -0.3:
        # if abs(tmp_particle.r - particle.r) < 0.3 * constants.a_0:
        if abs(tmp_particle.r - particle.r) < 0.35 * constants.a_0:
            # print(abs(tmp_particle.r - particle.r) / constants.a_0)
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
        u = np.random.uniform(0, 2 * np.pi)
        # v = np.random.uniform(0, np.pi)

        # cosu = np.random.uniform(-1, 1)
        cosv = np.random.uniform(-1, 1)
        tmp_z = np.random.uniform(0, 1)

        v = np.arccos(cosv)

        mul = tmp_z ** (1 / 3)
        r1 = radii[0] * mul
        r2 = radii[1] * mul
        r3 = radii[2] * mul

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


def get_and_plot_density(particles, radii, rot_axis, num=100):
    particles_ = copy.deepcopy(particles)

    for i in range(len(particles)):
        particles_[i].r.x *= 1 / radii[0]
        particles_[i].r.y *= 1 / radii[1]
        particles_[i].r.z *= 1 / radii[2]

    lengths = sorted([abs(particle.r) for particle in particles_])
    density = []
    l = []
    av_omega_list = []
    av_random_v_m_list = []
    particles = sorted(particles, key=lambda particle: abs(particle.r))

    for i in range(len(lengths) // num):
        first_ind = i * num
        last_ind = i * num + num - 1

        density.append(num * particles_[0].mass / (4 / 3 * np.pi * (lengths[last_ind] ** 3 - lengths[first_ind] ** 3)))
        # l.append(lengths[i*num] + 0.5*(lengths[i*num + num-1] - lengths[i*num]))
        l.append(lengths[first_ind])

        av_omega = 0
        for p in particles[first_ind:last_ind + 1]:
            av_omega += get_omega(p, rot_axis)
        av_omega /= len(particles[first_ind:last_ind + 1])
        av_omega_list.append(av_omega)

        av_rand_v = 0
        for p in particles[first_ind:last_ind + 1]:
            tmp_v = get_vel(p, rot_axis, av_omega)
            av_rand_v += abs(p.v - tmp_v)
        av_rand_v /= len(particles[first_ind:last_ind + 1])
        av_random_v_m_list.append(av_rand_v)

    first_ind = len(lengths) // num * num
    last_ind = -1

    density.append(num * particles_[0].mass / (4 / 3 * np.pi * (lengths[last_ind] ** 3 - lengths[first_ind] ** 3)))
    # l.append(lengths[i*num] + 0.5*(lengths[i*num + num-1] - lengths[i*num]))
    l.append(lengths[first_ind])

    av_omega = 0
    for p in particles[first_ind:last_ind]:
        av_omega += get_omega(p, rot_axis)
    av_omega /= len(particles[first_ind:last_ind])
    av_omega_list.append(av_omega)

    av_rand_v = 0
    for p in particles[first_ind:last_ind]:
        tmp_v = get_vel(p, rot_axis, av_omega)
        av_rand_v += abs(p.v - tmp_v)
    av_rand_v /= len(particles[first_ind:last_ind])
    av_random_v_m_list.append(av_rand_v)

    l.append(lengths[-1])

    # ig, ax1 = plt.subplots(figsize=(4, 4))
    #
    # ax1.scatter(x=l, y=density, marker='o', c='r', edgecolor='b')
    # ax1.set_title('Scatter: $x$ versus $y$')
    # ax1.set_xlabel('$x$')
    # ax1.set_ylabel('$y$')
    # plt.show()

    # io_xyz.write(particles_, save_path, -999)
    return l, av_random_v_m_list, av_omega_list, num


def create_sphere2(radii, l, num_sumples=2000, mass=1):
    attempts = 100
    # particles = [Particle(r=Vector3d(), mass=10)]
    particles = []
    num = 100
    for i, r in enumerate(l[1:]):
        # print('o', i)
        cur_num_particles = 0
        while True:
            # print(cur_num_particles)
            u = np.random.uniform(0, 2 * np.pi)
            cosv = np.random.uniform(-1, 1)
            tmp_z = np.random.uniform(0, 1)
            v = np.arccos(cosv)

            mul = tmp_z ** (1 / 3)

            mul = l[i] + mul * (r - l[i])

            r1 = radii[0] * mul
            r2 = radii[1] * mul
            r3 = radii[2] * mul

            tmp = [r1 * np.cos(u) * np.sin(v), r2 * np.sin(u) * np.sin(v), r3 * np.cos(v)]
            if criterion(tmp, particles):
                [x, y, z] = [tmp[0], tmp[1], tmp[2]]
                particles.append(Particle(r=Vector3d(x, y, z), mass=mass))

                cur_num_particles += 1

            if cur_num_particles == num // 1:  ##!!!!
                break
    return particles


def create_sphere3(radii, num_sumples=2000, mass=1):
    cur_num_particles = 0
    particles = []
    while True:
        u = np.random.uniform(0, 2 * np.pi)
        cosv = np.random.uniform(-1, 1)
        tmp_z = np.random.uniform(0, 1)
        v = np.arccos(cosv)

        mul = tmp_z ** (1 / 3)

        r1 = radii[0] * mul
        r2 = radii[1] * mul
        r3 = radii[2] * mul

        tmp = [r1 * np.cos(u) * np.sin(v), r2 * np.sin(u) * np.sin(v), r3 * np.cos(v)]
        if criterion(tmp, particles):
            [x, y, z] = [tmp[0], tmp[1], tmp[2]]
            particles.append(Particle(r=Vector3d(x, y, z), mass=mass))

            cur_num_particles += 1

        if cur_num_particles == num_sumples:  ##!!!!
            break
    return particles


def create_sphere4(radii, info, rot_axis, num_sumples=2000, mass=1):
    particles = []
    l, av_random_v_m_list, av_omega_list, num = info
    for i, el in enumerate(l[:-1]):
        cur_num_particles = 0
        r_min = el
        r_max = l[i + 1]
        av_random_v_m = av_random_v_m_list[i]
        av_omega = av_omega_list[i]
        while True:
            u = np.random.uniform(0, 2 * np.pi)
            cosv = np.random.uniform(-1, 1)
            tmp_z = np.random.uniform(0, 1)
            v = np.arccos(cosv)

            mul = tmp_z ** (1 / 3)
            mul = (r_max - r_min) * mul + r_min

            r1 = radii[0] * mul
            r2 = radii[1] * mul
            r3 = radii[2] * mul

            tmp = [r1 * np.cos(u) * np.sin(v), r2 * np.sin(u) * np.sin(v), r3 * np.cos(v)]

            if criterion(tmp, particles):
                [x, y, z] = [tmp[0], tmp[1], tmp[2]]
                tmp_par = Particle(r=Vector3d(x, y, z), mass=mass)
                tmp_par.v = get_vel(tmp_par, rot_axis, av_omega)

                tmp_e = Vector3d(np.random.sample(), np.random.sample(), np.random.sample())
                tmp_e *= (1 / abs(tmp_e))

                tmp_par.v += tmp_e * av_random_v_m

                particles.append(tmp_par)
                cur_num_particles += 1

            if (cur_num_particles == num) or ((cur_num_particles == num_sumples - (len(l)-2)*num) and (i == len(l) - 2)):  ##!!!!
                break

    mean_vel = Vector3d()
    for particle in particles:
        mean_vel += particle.v
    mean_vel *= (1 / len(particles))
    for particle in particles:
        particle.v -= mean_vel

    return particles


input_path = os.path.dirname(os.path.realpath(__file__)) + '/clusters/clust1804.xyz'
# input_path = '../3d/test/results_-999.xyz'
save_path = os.path.dirname(os.path.realpath(__file__)) + '/a_1957_thick_no_cheat_365'

if not os.path.exists(save_path):
    os.makedirs(save_path)

particles = io_xyz.read(input_path, mode='Nick')
a = constants.a_0

potential_energy = 0
for i, particle1 in enumerate(particles):
    for particle2 in particles[i + 1:]:
        potential_energy += forces_a.potential(particle1, particle2, a)
print('p', potential_energy)

print('Number of particles before: ', len(particles))

(center, radii, rotation) = getMinVolEllipse(
    np.array([[particle.r.x, particle.r.y, particle.r.z] for particle in particles], dtype='float32'))

mean_vel = Vector3d()
for particle in particles:
    mean_vel += particle.v
mean_vel = mean_vel * (1 / len(particles))

for particle in particles:
    particle.v -= mean_vel

kinetic_energy = 0
for particle in particles:
    kinetic_energy += particle.mass * abs(particle.v) ** 2 / 2

for particle in particles:
    particle.r -= Vector3d(center[0], center[1], center[2])
    tmp = np.dot(rotation, [particle.r.x, particle.r.y, particle.r.z])
    particle.r = Vector3d(tmp[0], tmp[1], tmp[2])

    tmp = np.dot(rotation, [particle.v.x, particle.v.y, particle.v.z])
    particle.v = Vector3d(tmp[0], tmp[1], tmp[2])

# find rot axis
tmp = np.zeros(3)
tmp[radii.argmin(axis=0)] = 1
rot_axis = Vector3d(tmp[0], tmp[1], tmp[2])

info = get_and_plot_density(particles, radii, rot_axis, 250)
print(len(info[1]))
io_xyz.write(particles, save_path, -1)
#################################################################################

particle_params_ratio = 1
particles = create_sphere4(radii, info, rot_axis, len(particles) // particle_params_ratio, mass=particle_params_ratio)

mi_r = 1000
for i, p1 in enumerate(particles):
    for p2 in particles[i + 1:]:
        if mi_r > abs(p2.r - p1.r):
            mi_r = abs(p2.r - p1.r)
print(mi_r)
io_xyz.write(particles, save_path, -2)
print(len(particles))

(center1, radii1, rotation1) = getMinVolEllipse(
    np.array([[particle.r.x, particle.r.y, particle.r.z] for particle in particles], dtype='float32'))
print(center1, radii1, rotation1)

potential_energy_new = 0
for i, particle1 in enumerate(particles):
    for particle2 in particles[i + 1:]:
        potential_energy_new += forces_a.potential(particle1, particle2, a)
print('pe', potential_energy - potential_energy_new, potential_energy / potential_energy_new)

for i in range(1000):
    # a += -0.001/(i+1)*np.sign(potential_energy - potential_energy_new)
    a += -0.001/((potential_energy - potential_energy_new)/potential_energy_new)

    potential_energy_new = 0
    for i, particle1 in enumerate(particles):
        for particle2 in particles[i + 1:]:
            potential_energy_new += forces_a.potential(particle1, particle2, a)
    print('pe', potential_energy - potential_energy_new, potential_energy / potential_energy_new)

    if abs(1 - potential_energy / potential_energy_new) < 0.1:
        break

print("a", a)

kinetic_energy_new = 0
for particle in particles:
    kinetic_energy_new += particle.mass * abs(particle.v) ** 2 / 2
print("ke", kinetic_energy - kinetic_energy_new, kinetic_energy / kinetic_energy_new)

#####
for step in range(constants.steps_number):
    print('s', step)
    if step % 10 == 0:
        io_xyz.write(particles, save_path, step)

    for particle in particles:
        particle.force = Vector3d()

    for i, p1 in enumerate(particles):
        for p2 in particles[i + 1:]:
            tmp_force = forces_a.base_force(p1, p2, a)
            p1.force += tmp_force
            p2.force += -tmp_force

    for particle in particles:
        particle.v += particle.force * constants.dt
        particle.r += particle.v * constants.dt
#
