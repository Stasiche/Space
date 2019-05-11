from vector3d import Vector3d


class Particle:
    def __init__(self, r, name='noname', v=Vector3d(0, 0), mass=1, color=(1, 1, 1)):
        self.r = r
        self.v = v
        self.name = name
        self.mass = mass
        self.force = Vector3d(0, 0)
        self.color = color

    def add_r(self, other):
        self.r += other.r

    def add_v(self, other):
        self.v += other.v

    def set_v(self, new_v):
        self.v = new_v

    def __str__(self):
        return 'r: ' + str(self.r) + '\nv: ' + str(self.v)