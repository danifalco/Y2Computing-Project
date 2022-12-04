import matplotlib.pyplot as plt

from raytracer import *
from copy import deepcopy

ray_list = col_generator()

ray_list_mess = deepcopy(ray_list)
# ray_list = [Ray(np.array([i, 0, 0]), np.array([0, 0, 1])) for i in np.arange(-5, 5.5, 0.5)]
plano_conv = BiConvex(101, 0.02, -0.02, 1.5, 5, 1)
# out_pln = OutputPlane(1500, 500)
plano_conv.plot_xz(ray_list_mess)


# for ray in ray_list_mess:
#     points = ray.vertices()
#     plt.plot([i[2] for i in points],
#              [i[0] for i in points], '-x', color='k')
plt.title('hola')
plt.show()

# curvature = 0.02
#
# spherical = SphericalRefraction(100, -curvature, 1.5, 1, 5)
# # spherical2 = SphericalRefraction(101, -curvature, 1.5, 1, 5)
# for ray in ray_list:
#     spherical.propagate_ray(ray)
#     # spherical2.propagate_ray(ray)
#     out_pln.propagate_ray(ray)
#     points = ray.vertices()
#     plt.plot([i[2] for i in points],
#              [i[0] for i in points], '-x', color='k')

    # xz_points = [[i[2], i[0]] for i in points]
    # plt.scatter(xz_points[0][0], xz_points[0][1], color='g')
    # plt.scatter(xz_points[1][0], xz_points[1][1], color='b')
    # plt.scatter(xz_points[2][0], xz_points[2][1], color='r')
# plt.title(curvature)
# plt.xlim(95, 10)
# plt.show()

######

z, rms = rms_plotter(ray_list_mess, np.arange(127.5, 135, 0.05))
print(f"Min RMS is in position {z[rms.index(min(rms))]}")
plt.plot(z, rms, '-x')
plt.show()
