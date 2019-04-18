# particles - [0rx,1ry,2rz,3vx,4vy,5vz,6mass,7fx,8fy,9fz]
import numpy as np
import copy

a = 0.01
theta = 0.0
number_particles = 2
number_steps = 1000
dt = 0.00001

k = 25
l = 25
calculate_area = a * k * 4

saved_particles = 0
min_mac = 100000000

save_path = '/Users/stas/PycharmProjects/Space/fast/'


def make_hex(k, l):
    result = np.array([])
    a_sqrt3_2 = a * np.sqrt(3) * 0.5
    for i in range(k):
        for j in range(i + l):
            particle = np.array([])
            particle = np.append(particle, i*a_sqrt3_2)
            particle = np.append(particle, a*j - a * 0.5 * i)
            particle = np.append(particle, 0)

            for iv in range(3):
                particle = np.append(particle, 0)

            particle = np.append(particle, 1)

            for fi in range(3):
                particle = np.append(particle, 0)

            if len(result) == 0:
                # result = copy.deepcopy(particle)
                result = particle
            else:
                result = np.vstack((result, particle))

    for i in range(k, 2 * k - 1):
        for j in range(l + k - 2 - (i-k)):
            particle = np.array([])
            particle = np.append(particle, i*a_sqrt3_2)
            particle = np.append(particle, a*j + a * 0.5 * (i - 2 * k + 2))
            particle = np.append(particle, 0)

            for iv in range(3):
                particle = np.append(particle, 0)

            particle = np.append(particle, 1)

            for fi in range(3):
                particle = np.append(particle, 0)

            result = np.vstack((result, particle))

    return result


def force(particle1, particle2):
    l = 0
    # p1 = copy.deepcopy(particle1)
    # p2 = copy.deepcopy(particle2)
    p1 = particle1
    p2 = particle2
    e = np.zeros(3)
    result = np.zeros(3)

    for i in range(3):
        e[i] = p2[i] - p1[i]

    for i in range(3):
        l += e[i] ** 2
    l = np.sqrt(l)

    for i in range(3):
        e[i] /= l

    for i in range(3):
        result[i] = e[i] * (1 / a**2) * (a/l) ** 2 * (1 - (a/l)**4)

    return result


def save_xyz_r(name, particles):
    with open(name + '.xyz', 'w') as outfile:
        outfile.write(str(len(particles)) + '\n\n')
        for i, particle in enumerate(particles):
            outfile.write(str(particle[0]) + ' ' +
                          str(particle[1]) + ' ' +
                          str(particle[2]) + ' ' + '\n')


particles = make_hex(k, l)

for step in range(number_steps):
    # clear forces
    for i in range(len(particles)):
        for fi in range(3):
            particles[i][7 + fi] = 0

    # calculate forces
    for i in range(len(particles)):
        # print(i)
        for j in range(i+1, len(particles)):
            brute_force = force(particles[i], particles[j])
            for fi in range(3):
                particles[i][7 + fi] += brute_force[fi]
                particles[j][7 + fi] += -brute_force[fi]


     # move particles
    for i in range(len(particles)):
        for j in range(3):
            particles[i][3 + j] += particles[i][7 + j] * dt

        for j in range(3):
            particles[i][j] += particles[i][3 + j] * dt

    # save particles
    save_xyz_r(save_path + str(step), particles)
    if step % 1 == 0:
        print(step)
