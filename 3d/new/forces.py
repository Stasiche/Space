import constants
from vector3d import Vector3d


def base_force(particle1, particle2):
    r_abs = abs(particle1.r - particle2.r)
    e = (particle1.r - particle2.r) * (1/r_abs)

    f_tmp = particle2.mass * particle1.mass * constants.gamma_0 * 1/(r_abs ** 2)

    if (r_abs ** 2) < constants.a_cut_sqr_0:
        dvr = (particle2.v - particle1.v).dot(particle2.r - particle1.r)
        t4 = (constants.a_0 / r_abs) ** 4
        Q = constants.beta * dvr * t4 / (r_abs ** 2)
        # f_tmp *= (1 - t4 + Q)
        f_tmp *= (1 - t4)
    # return -e * (gamma * 1/(constants.a ** 2)*(constants.a/r_abs)**2*(1-(constants.a/r_abs)**4))
    return -e * f_tmp


def potential(particle1, particle2):
    r_abs = abs(particle1.r - particle2.r)

    p_tmp = particle2.mass * particle1.mass * constants.gamma_0 * 1/(r_abs ** 2)

    if (r_abs ** 2) < constants.a_cut_sqr_0:
        p_tmp -= particle1.mass*particle2.mass*constants.gamma_0*(constants.a_0)**4/(5*r_abs**5)
    return p_tmp

