from particle import Particle
from vector2d import Vector2d


def write(particles, save_path, step):
    with open(save_path + '/results_' + str(step) + '.xyz', 'w') as outfile:
        outfile.write(str(len(particles)) + '\n\n')
        for i, particle in enumerate(particles):
            outfile.write(str(particle.r) + ' ' +
                          str(particle.color[0]) + ' ' +
                          str(particle.color[1]) + ' ' +
                          str(particle.color[2]) + ' ' + '\n')


def read(read_path):
    particles = []
    with open(read_path, 'r') as inputfile:
        s_tmp = inputfile.read()
        for el in s_tmp.split('\n')[2:]:
            if el != '':
                rx, ry, cr, cg, cb = el.split()
                particles.append(Particle(r=Vector2d(float(rx), float(ry)), name=len(particles)))
    return particles