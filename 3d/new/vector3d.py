from math import sqrt


class Vector3d:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def zeros(self):
        self.x = 0
        self.y = 0
        self.z = 0

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def __iadd__(self, other):
        return Vector3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __neg__(self):
        return self.__mul__(-1)

    def __add__(self, other):
        return Vector3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vector3d(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        return Vector3d(other*self.x, other*self.y, other*self.z)

    def __str__(self):
        return str(self.x) + ' ' + str(self.y) + ' ' + str(self.z)

    def __abs__(self):
        return sqrt(self.x**2 + self.y**2 + self.z**2)
