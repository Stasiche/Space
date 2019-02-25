from vector2d import Vector2d


class Node:
    # ROOT = 0
    # INTERNAL = 1
    # EXTERNAL = 2

    def __init__(self, rect, parent=None, node_type=None, name=None):
        self.parent = parent
        self.children = [None, None, None, None]

        if parent is None:
            self.depth = 0
            self.name = 'n'
        else:
            self.depth = parent.depth + 1
            self.name = self.parent.name + '_' + name

        self.rect = rect

        # if self.parent is None:
        #     self.type = 0
        # else:
        #     self.type = node_type
        self.type = node_type
        self.particles = []

        self.mass = None
        self.cm = None                  # non-normed!
        self.mv = None                  # non-normed!

    def add_particle(self, particle):
        self.particles.append(particle)
        self.update_mass_cm(particle)

    def update_mass_cm(self, particle):
        if self.mass is None:
            self.mass = particle.mass
        else:
            self.mass += particle.mass

        if self.cm is None:
            self.cm = particle.r * particle.mass
        else:
            self.cm += particle.r * particle.mass

        if self.mv is None:
            self.mv = particle.v
        else:
            self.mv += particle.v

    def get_cm(self):
        return self.cm * (1/self.mass)

    def get_mv(self):
        return self.mv * (1/len(self.particles))

    def check_particle_in_node(self, particle):
        x0, y0, x1, y1 = self.rect
        if (particle.r.x >= x0) and (particle.r.y >= y0) and (particle.r.x < x1) and (particle.r.y < y1):
            return True
        else:
            return False

    def subdivide(self):
        x0, y0, x1, y1 = self.rect
        h = (x1 - x0) / 2
        rects = list()
        rects.append((x0, y0, x0 + h, y0 + h))
        rects.append((x0, y0 + h, x0 + h, y1))
        rects.append((x0 + h, y0 + h, x1, y1))
        rects.append((x0 + h, y0, x1, y0 + h))

        for i, rect in enumerate(rects):
            self.children[i] = Node(parent=self, rect=rect, node_type=2, name=str(i))

        self.type = 1

