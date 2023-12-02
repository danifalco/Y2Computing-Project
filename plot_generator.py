import matplotlib.pyplot as plt
from raytracer import *
from copy import deepcopy

# %% Plot 1 - From Task 9

ray_list = [Ray(np.array([i, 0, 0]), np.array([0, 0, 1])) for i in np.arange(-5, 5.5, 0.5)]

ref_surface = SphericalRefraction(100, 0.03, 1, 1.5, 10)
out_pln = OutputPlane(250, 5)
for ray in ray_list:
    ref_surface.propagate_ray(ray)
    out_pln.propagate_ray(ray)
    points = ray.vertices()
    plt.plot([i[2] for i in points],
             [i[0] for i in points], color='#1f77b4')

plt.title("Ray Trajectory Through Spherical Surface")
plt.xlabel("Z-Axis (mm)")
plt.ylabel("X-Axis (mm)")
plt.grid()
plt.savefig("Plots\\Plot 1", dpi=700)
plt.show()

# %% Plot 2 - Task 12
ray_list = col_generator()
ref_surface = SphericalRefraction(100, 0.03, 1, 1.5, 10)
out_pln = OutputPlane(250, 5)
for ray in ray_list:
    ref_surface.propagate_ray(ray)
    out_pln.propagate_ray(ray)
    points = ray.vertices()
    plt.plot([i[2] for i in points],
             [i[0] for i in points], color='#1f77b4')

plt.title("Large Diameter Beam Refraction")
plt.xlabel("Z-Axis (mm)")
plt.ylabel("X-Axis (mm)")
plt.grid()
plt.savefig("Plots\\Plot 2", dpi=700)
plt.show()

# %% Plot 3 - Task 13

ray_list = col_generator()
ref_surface = SphericalRefraction(100, 0.03, 1, 1.5, 10)
out_pln = OutputPlane(200, 5)
plt.figure(figsize=(8, 8))
for ray in ray_list:
    ref_surface.propagate_ray(ray)
    out_pln.propagate_ray(ray)
    points = ray.p()
    plt.scatter(points[0], points[1], s=7,color='#1f77b4')

plt.title("Spot Diagram (z = 200mm)")
plt.xlabel("X-Axis (mm)")
plt.ylabel("Y-Axis (mm)")
plt.savefig("Plots\\Plot 3", dpi=700)
plt.grid(True)
plt.show()

# %% Plot 4
ray_list = col_generator()

plano_convex_1 = PlanoConvex(100, 0, -0.02, 1.5, 7)  # Flat side first
plano_convex_2 = PlanoConvex(100, 0.02, 0, 1.5, 7)  # Flat side last
biconvex = BiConvex(100, 0.02, -0.02, 1.5, 7)
out_pln = OutputPlane(180, 5)

ray_list_copy = deepcopy(ray_list)
plano_convex_1.plot_xz(ray_list_copy, out_pln)
plt.title("Plano-Convex Refraction - Flat Side First")
plt.xlabel("Z-Axis (mm)")
plt.ylabel("X-Axis (mm)")
plt.grid()
plt.savefig("Plots\\Plot 4a", dpi=700)
plt.show()

ray_list_copy = deepcopy(ray_list)
plano_convex_2.plot_xz(ray_list_copy, out_pln)
plt.title("Plano-Convex Refraction - Curved Side First")
plt.xlabel("Z-Axis (mm)")
plt.ylabel("X-Axis (mm)")
plt.grid()
plt.savefig("Plots\\Plot 4b", dpi=700)
plt.show()

ray_list_copy = deepcopy(ray_list)
out_pln = OutputPlane(200, 5)
biconvex.plot_xz(ray_list_copy, out_pln)
plt.title("Biconvex Lens Refraction")
plt.xlabel("Z-Axis (mm)")
plt.ylabel("X-Axis (mm)")
plt.grid()
plt.savefig("Plots\\Plot 4c", dpi=700)
plt.show()


# %% Plot 5

plano_convex_1 = PlanoConvex(100, 0, -0.02, 1.5, 22)  # Flat side first
plano_convex_2 = PlanoConvex(100, 0.02, 0, 1.5, 22)  # Flat side last
biconvex = BiConvex(100, 0.02, -0.02, 1.5, 22)

out_pln = OutputPlane(154.89, 5)
for i in np.arange(0.1, 20, 0.3):
    ray_list = col_generator(rad=i, dist_pts=0.3)
    for ray in ray_list:
        plano_convex_1._propagate(ray)
    plt.scatter(i, rms(ray_list), color='r', s=5)
plt.scatter(0, 0, color='r', s=5, label="Plano-Convex - Curved Side First")

out_pln = OutputPlane(177.41, 5)
for i in np.arange(0.1, 20, 0.3):
    ray_list = col_generator(rad=i, dist_pts=0.3)
    for ray in ray_list:
        plano_convex_2._propagate(ray)
    plt.scatter(i, rms(ray_list), color='b', s=5)
plt.scatter(0, 0, color='b', s=5, label="Plano-Convex -  Flat Side First")

out_pln = OutputPlane(139.25, 5)
for i in np.arange(0.1, 20, 0.3):
    ray_list = col_generator(rad=i, dist_pts=0.3)
    for ray in ray_list:
        biconvex._propagate(ray)
    plt.scatter(i, rms(ray_list), color='g', s=5)
plt.scatter(0, 0, color='g', s=5, label="Biconvex")

plt.xlabel("Beam radius (mm)")
plt.ylabel("RMS Radius (mm)")
plt.title("Performance of Different Lenses")
plt.legend(loc='lower right')
plt.grid()
plt.savefig("Plots\\Plot 5", dpi=700)
plt.show()
