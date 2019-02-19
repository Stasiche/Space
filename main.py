from particle import Particle
from vector2d import Vector2d

c = 1e0
l0 = 10
t = 1e5
dt = 1/1e0

grav_coef = 1/1e0

par_list = []
par_list.append(Particle(r=Vector2d(10, 0), v=Vector2d(0, 0)))
par_list.append(Particle(r=Vector2d(20, 0.01), v=Vector2d(0, 0)))


# while True:
while t > 0:
    # l_v = par_list[0].r - par_list[1].r
    # l = abs(l_v.__abs__())
    # inv_l = 1/l
    # e_v = -l_v*inv_l

    l_v = par_list[0].r
    l = abs(l_v.__abs__())
    inv_l = 1/l
    e_v = -l_v*inv_l
    # print(e_v.__abs__())

    print(l)
    # f_upr = l*dt*c
    f_upr = 0
    f_grav_0 = e_v * grav_coef * inv_l**2 * dt
    # print(f_upr)

    # f_0 = f_upr + f_grav_0
    f_0 = f_grav_0
    print(f_0)
    par_list[0].v += f_0
    print(par_list[0].v)
    # par_list[1].v += -e_v*f_upr

    par_list[0].r += par_list[0].v*dt
    # par_list[1].r += par_list[1].v*dt

    t -= 1
    print('_')