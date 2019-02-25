import pylab
import constants


def draw_rect(axes, x0, y0, x1, y1):
    rect_coord = (x0, y0)
    rect_width = x1 - x0
    rect_height = y1 - y0

    rect = pylab.Rectangle(rect_coord, rect_width, rect_height, color='#000000', fill=False)
    axes.add_patch(rect)


def plot_nodes(axes, node):
    draw_rect(axes, *node.rect)
    for child in node.children:
        if child is not None:
            plot_nodes(axes, child)


def plot_particles(axes, particles):
    x = [particle.r.x for particle in particles]
    y = [particle.r.y for particle in particles]
    # axes.scatter(x=x, y=y, marker='o', c='r', edgecolor='b')
    # axes.set_title('Scatter: $x$ versus $y$')
    # axes.set_xlabel('$x$')
    # axes.set_ylabel('$y$')

    # pylab.scatter(x, y, c=color, s=size, alpha=alpha, edgecolors='none')
    pylab.scatter(x, y, s=constants.a / 4)


def plot_main(particles, bh_tree):
    x0, y0, x1, y1 = bh_tree.root.rect
    pylab.xlim(x0, x1)
    pylab.ylim(y0, y1)

    # Получим текущие оси
    axes = pylab.gca()
    axes.set_aspect("equal")

    plot_nodes(axes, bh_tree.root)
    plot_particles(axes, particles)

    pylab.show()
    # pylab.ion()

