import constants


def mac(particle, node):
    x0, y0, x1, y1 = node.rect
    px = particle.r.x
    py = particle.r.y
    cm = node.get_cm()
    # return abs(x1-x0) / (abs(particle.r - cm) + constants.epsilon) < constants.theta
    return abs(x1-x0) / (min([abs(x0 - px), abs(x1 - px),
                              abs(y0 - py), abs(y1 - py)]) + constants.epsilon) < constants.theta
