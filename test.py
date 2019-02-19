from node import Node
from particle import Particle
from vector2d import Vector2d
from bh_tree import BHtree
from profiler import Profiler
import random
import forces
import constants

import networkx as nx
import matplotlib.pyplot as plt

def make_univer(n_particles=int(10e6)):
    max_v = 1
    min_v = -1
    particles = [Particle(r=Vector2d(random.uniform(-box_size // 2, box_size // 2),
                                             random.uniform(-box_size // 2, box_size // 2)),
                                  v=Vector2d(random.uniform(min_v, max_v),
                                             random.uniform(min_v, max_v)))]
    # global box_size

    i = 1
    while True:
        # if i % (n_particles//100) == 0:
        #     print('box:', i // (n_particles//100) + 1)
        tmp_par = Particle(r=Vector2d(random.uniform(-box_size // 2, box_size // 2),
                                             random.uniform(-box_size // 2, box_size // 2)),
                                  v=Vector2d(random.uniform(min_v, max_v),
                                             random.uniform(min_v, max_v)))
        flag = 1
        # for particle1 in particles:
        #     tmp_e = particle1.r - tmp_par.r
        #     if abs(tmp_e) == 0:
        #         continue
        #     tmp_f = forces.base_force(particle1, tmp_par)
        #     if (tmp_f.x*tmp_e.x + tmp_f.y*tmp_e.y) < 0:
        #         flag = 0

        if flag == 1:
            particles.append(tmp_par)
            i += 1

        if i >= n_particles:
            break

        # print(i)
    return particles


with Profiler() as p:
    box_size = constants.box_size
    npar = constants.npar
    particles = make_univer(npar)

    for t in range(constants.steps_number):
        with open('results_' + str(t) + '.xyz', 'w') as outfile:
            outfile.write(str(len(particles)) + '\n\n')
            for i, particle in enumerate(particles):
                # if i % (npar // 100) == 0:
                #     print('write:', i // (npar // 100) + 1)
                outfile.write(str(particle.r)+' 250 \n')

        method = constants.method
        if method == 'bh':
            # Barnes-Hut #########################
            n1 = Node([-box_size / 2, -box_size / 2, box_size / 2, box_size / 2], node_type=2)
            bh = BHtree(n1)
            for i, particle in enumerate(particles):
                # if i % (npar // 100) == 0:
                #     print('tree:', i // (npar // 100) + 1)
                bh.insert_particle(particle, bh.root)

            # add_list = bh.gen_add_list(bh.root, [])
            # G = nx.Graph()
            # G.add_edges_from(add_list)
            # nx.draw_networkx(G)
            # plt.show()
            # bh.count_external_nodes(bh.root)
            # print(bh.external_nodes_num)

            for i, particle in enumerate(particles):
                # if i % (npar // 100) == 0:
                #     print('force:', i // (npar // 100) + 1)
                particle.force = bh.calculate_force(particle, bh.root)
            #########################################
        elif method == 'brute':
            # Brute ###################
            for i, particle1 in enumerate(particles):
                particle1.force = Vector2d(0, 0)

            for i, particle1 in enumerate(particles):
                for particle2 in particles[i+1:]:
                    tmp_force = forces.base_force(particle1, particle2)
                    particle1.force += tmp_force
                    particle2.force -= tmp_force
                #
                # if i % (npar // 100) == 0:
                #     print('brute:', i // (npar // 100) + 1)
            #############################

        for i, particle in enumerate(particles):
            # if i % (npar // 100) == 0:
            #     print('move:', i // (npar // 100) + 1)

            # particle.v -= particle.r * constants.dt
            particle.v += -particle.force * constants.dt
            particle.r += particle.v * constants.dt

        else:
            mean_error = 0
            n1 = Node([-box_size / 2, -box_size / 2, box_size / 2, box_size / 2], node_type=2)
            bh = BHtree(n1)
            for i, particle in enumerate(particles):
                bh.insert_particle(particle, bh.root)

            for i, particle in enumerate(particles):
                bh_force = bh.calculate_force(particle, bh.root)
                # particle.force = bh_force

                brute_force = Vector2d(0, 0)
                for particle2 in particles:
                    if particle != particle2:
                        brute_force += forces.base_force(particle, particle2)

                particle.force = bh_force

                print(abs(bh_force - brute_force))
                mean_error += abs(bh_force - brute_force)
            mean_error /= len(particles)
            #########################################

        print(t, 'done!', '\nmean error: ' + str(mean_error))
    print('done!!')
