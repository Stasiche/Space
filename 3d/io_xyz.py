from particle import Particle
from vector3d import Vector3d


def write(particles, save_path, step):
    with open(save_path + '/results_' + str(step) + '.xyz', 'w') as outfile:
        outfile.write(str(len(particles)) + '\n\n')
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.color[0]) + ' ' +
                          str(particle.color[1]) + ' ' +
                          str(particle.color[2]) + ' ' + '\n')


def append_write(particles, save_path, step):
    particles += read(save_path + '/results_' + str(step) + '.xyz')
    with open(save_path + '/results_' + str(step) + '.xyz', 'w') as outfile:
        outfile.write(str(len(particles)) + '\n\n')
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.color[0]) + ' ' +
                          str(particle.color[1]) + ' ' +
                          str(particle.color[2]) + ' ' + '\n')


def read(read_path, mode=None):
    particles = []
    with open(read_path, 'r') as inputfile:
        s_tmp = inputfile.read()
        for el in s_tmp.split('\n')[2:]:
            if el != '':
                if mode is None:
                    rx, ry, rz, r, g, b = el.split()
                else:
                    с, rx, ry, rz, c = el.split()
                particles.append(Particle(r=Vector3d(float(rx), float(ry), float(rz)), name=len(particles)))
    return particles