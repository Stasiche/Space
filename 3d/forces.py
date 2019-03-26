import constants
from vector2d import Vector2d


def base_force(particle1, particle2):
    r_abs = abs(particle1.r - particle2.r)
    e = (particle1.r - particle2.r) * (1/r_abs)
    gamma = particle2.mass * particle1.mass

    return -e * (gamma * 1/(constants.a ** 2)*(constants.a/r_abs)**2*(1-(constants.a/r_abs)**4))
    # return e * (12*constants.D/constants.a*(1-(constants.a/r_abs)**6))
    # return e * (constants.G / (1-(constants.a/r_abs) ** 6))
    # return e * (constants.G * particle1.mass * particle2.mass / (r_abs ** 3))

    # return e * (12*constants.D/constants.a*(constants.a/r_abs)**7*(1-(constants.a/r_abs)**6))

