import numpy as np

from raytracer import *
import matplotlib.pyplot as plt
import copy
########################################################################################################################
# snell_refraction(np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0]), np.array([1, 0, 0]), 1.0, 1.5)
#
# ray1 = Ray(np.array([1, 0, -10]), np.array([0, 0, 2]))
#
# ray1.append(np.array([1, 2, 3]), np.array([1, 2, 3]))
#
# sphere = SphericalRefraction(0, 0.2, 1, 1.5, 0)

# print(sphere._intercept(ray1))
# sphere.propagate_ray(ray1)
# print(ray1.is_terminated())
# out_pln = OutputPlane(0, 5)
# print(out_pln.intercept_out(ray1))


########################################################################################################################
#
# print(sphere.intercept(ray1))
# print(type(thing.get_position()))

########



# raylist = [Ray(np.array([i / 150, 0, -100]), np.array([0, 0, 1])) for i in range(-15, 16)]




# for ray in raylist:
#     # ref_surface.intercept(ray)
#     ref_surface.propagate_ray(ray)
#     # out_pln.intercept(ray)
#     out_pln.propagate_ray(ray)
#     plt.plot([i[2] for i in ray.vertices()], [i[0] for i in ray.vertices()], color='k')
#
# plt.show()

# class ColimatedBeam:
#     def __init__(self, z0: float = 0, diam: int = 5, dpmm: int = 1, num_pts: int = 8):
#         self._z0 = z0

# def col_generator(z0: float = 0, rad: int = 5, dist_pts: float = 0.5, k_vec: Union[None, np.ndarray] = None) -> list:
#     if k_vec is None:  # Properly dealing with mutable objects passed as kwargs
#         k = np.array([0, 0, 1])
#     else:
#         k = k_vec
#     ray_lst = [Ray(np.array([0, 0, z0]), k)]
#
#     pts_incr = 6  # 6 is the best approximation to "uniformly distributed". This is impossible in reality
#     num_of_shells = int(rad / dist_pts)
#     for i in range(num_of_shells):
#         num_pts = i * pts_incr  # The number of rays in the "shell"
#         radius = i * dist_pts
#         ray_lst += (Ray(pos=np.array([radius * np.cos(2 * np.pi * dot / num_pts),
#                                       radius * np.sin(2 * np.pi * dot / num_pts),
#                                       z0]),
#                         direc=np.array([0, 0, 1]))
#                     for dot in range(num_pts))
#         # x_vals += (radius * np.cos(2 * np.pi * dot / num_pts) for dot in range(num_pts))
#         # y_vals += (radius * np.sin(2 * np.pi * dot / num_pts) for dot in range(num_pts))
#
#     # return x_vals, y_vals
#     return ray_lst



list_of_rays = col_generator()

# alst = [(ray.p()[0], ray.p()[1]) for ray in list_of_rays]
#
# plt.figure(figsize=(5, 5))
# for point in alst:
#     plt.scatter(point[0], point[1], color='r')

curvature = -0.02

ref_surface = SphericalRefraction(90, curvature, 1, 1.5, 3)
flat_surface = SphericalRefraction(100, 0, 1, 1, 3)
# out_pln = OutputPlane(150, 200)

for ray in list_of_rays:
    ref_surface.propagate_ray(ray)
    flat_surface.propagate_ray(ray)

    points = ray.vertices()
    plt.plot([i[2] for i in points],
             [i[0] for i in points], '-x', color='k')

plt.show()

ray_backup = copy.deepcopy(list_of_rays)

# for trial_z in np.arange(105, 185, 0.1):
#     trial_rays = copy.deepcopy(ray_backup)
#     out_pln = OutputPlane(trial_z, 200)
#     for ray in trial_rays:
#         out_pln.propagate_ray(ray)
#     plt.scatter(trial_z, rms(trial_rays), color='k')

z, rms = rms_plotter(ray_backup, np.arange(105, 185, 0.1))
plt.plot(z, rms, '-x')
plt.show()

    # out_pln.propagate_ray(ray)
    # plt.plot([i[2] for i in ray.vertices()],
    #             [i[1] for i in ray.vertices()], '-x', color='b')
    # plt.plot([i[2] if i[1] == 0 else None for i in ray.vertices()],
    #          [i[0] if i[1] == 0 else None for i in ray.vertices()])

# x = np.arange(100, 110, 0.1)
# plt.plot(x, np.sqrt(1 / curvature ** 2 - (x - 100) ** 2), color='r')
# position_list = [ray._position for ray in list_of_rays]
# xlist = [[]]


# for list_of_points in position_list:
#     plt.plot(list_of_points[2], list_of_points[0], color='b')
# plt.title(f'curvature = {curvature}')
# plt.show()
# "---------------"
# beam = col_generator()
#
# plano = SphericalRefraction(100, 0, 1, 1.5, 0)
# convex = SphericalRefraction(130, -0.2, 1.5, 1, 0)
# output_pln = OutputPlane(250)
#
# for ray in beam:
#     plano.propagate_ray(ray)
#     convex.propagate_ray(ray)
#     output_pln.propagate_ray(ray)
#     plt.plot(ray.p()[2], ray.p()[0])
#
# plt.show()



print("-----end of script-----")