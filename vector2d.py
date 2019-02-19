from math import sqrt


class Vector2d:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def zeros(self):
        self.x = 0
        self.y = 0

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def __iadd__(self, other):
        return Vector2d(self.x + other.x, self.y + other.y)

    def __neg__(self):
        return self.__mul__(-1)

    def __add__(self, other):
        return Vector2d(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Vector2d(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return Vector2d(other*self.x, other*self.y)

    def __str__(self):
        return str(self.x) + ' ' + str(self.y)

    def __abs__(self):
        return sqrt(self.x**2 + self.y**2)