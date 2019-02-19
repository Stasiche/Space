from vector2d import Vector2d
from particle import Particle

import forces
import constants

class BHtree:
    def __init__(self, root):
        self.root = root
        self.external_nodes_num = None

    def insert_particle(self, particle, node):
        if len(node.particles) == 0:
            node.add_particle(particle)

        elif node.type == 1:
            for child in node.children:
                if child.check_particle_in_node(particle):
                    self.insert_particle(particle, child)

            # node.update_mass_cm(particle)
            node.add_particle(particle)

        elif node.type == 2:
            node.add_particle(particle)
            node.subdivide()
            for internal_particle in node.particles:
                for child in node.children:
                    if child.check_particle_in_node(internal_particle):
                        self.insert_particle(internal_particle, child)

            # node.update_mass_cm(particle)

    def calculate_force(self, particle, node=None):
        if node is None:
            node = self.root

        force = Vector2d(0, 0)
        if len(node.particles) == 0:
            return force

        # if (node.type == 2) and (abs(node.particles[0].r - particle.r) != 0):
        if (node.type == 2) and (node.particles[0] != particle):
            tmp_f = forces.base_force(particle, node.particles[0])
            # print(particle.name, node.name, node.depth, len(node.particles))
            # if particle.name == 9:
            #     tmp_brute_f = Vector2d()
            #     for particle_brute in node.particles:
            #         tmp_brute_f += forces.base_force(particle, particle_brute)
            #
            #     len_error = (abs(tmp_f) - abs(tmp_brute_f)) / (abs(tmp_brute_f) + constants.epsilon)
            #     print(len_error, '|', tmp_f, '|', tmp_brute_f)
            force += tmp_f
            # force += forces.gravity(particle, node.particles[0])

        x0, y0, x1, y1 = node.rect
        cm = node.get_cm()
        if (node.type == 1) and (abs(x1-x0)/(abs(particle.r - cm) + constants.epsilon) < constants.theta):
            # print(particle.name, node.name)
            print(particle.name,  node.name, node.depth, len(node.particles))
            tmp_f = forces.base_force(particle, Particle(r=cm, mass=node.mass, name='aprox'))
            if particle.name == 9:
                tmp_brute_f = Vector2d()
                for particle_brute in node.particles:
                    tmp_brute_f += forces.base_force(particle, particle_brute)

                len_error = (abs(tmp_f) - abs(tmp_brute_f)) / (abs(tmp_brute_f) + constants.epsilon)
                print(len_error, '|', tmp_f, '|', abs(tmp_f), '|', tmp_brute_f, '|', abs(tmp_brute_f))
            force += tmp_f
            # force += forces.gravity(particle, Particle(r=cm, mass=node.mass))

        if (node.type == 1) and (abs(x1-x0)/(abs(particle.r - cm) + constants.epsilon) >= constants.theta):
            for child in node.children:
                force += self.calculate_force(particle, child)

        return force

    def gen_add_list(self, node, add_list):
        if node.type == 2:
            add_list.append((node.name, node.parent.name))
        else:
            tmp_add_list = []
            for child in node.children:
                tmp_add_list += self.gen_add_list(child, add_list)
            add_list += tmp_add_list
        return add_list

    def count_external_nodes(self, node):
        if self.external_nodes_num is None:
            self.external_nodes_num = 0

        if (node.type == 2) and (len(node.particles) != 0):
            self.external_nodes_num += 1
        elif node.type == 1:
            for child in node.children:
                self.count_external_nodes(child)
