'''
:Date: 17-11-2022
:Author: Daniel Weronski Falco, Blackett Laboratory

First module

Classes:


Functions:

'''

import numpy as np
from scipy.optimize import fmin_tnc
import matplotlib.pyplot as plt
import copy
import warnings

z_hat = np.array([0, 0, 1])


class Ray:
    """
    Ray class to represent a ray of light
    """

    def __init__(self, pos: np.ndarray, direc: np.ndarray) -> None:
        """
        Ray class constructor

        :param pos: vector space position in 3D: (x, y, z)
        :param direc: direction vector in 3D: (x, y, z)
        :return: None
        """
        self._position = np.array([pos])
        self._direction = np.array([direc])

    def __repr__(self):
        return f"Ray(pos=np.array([{self.p()[0]}, {self.p()[1]}, {self.p()[2]}]), direc=np.array([{self.k()[0]}, " \
               f"{self.k()[1]}, {self.k()[2]}]))"

    def __str__(self):
        return f"Ray object. Currently in position {self.p()} and with direction {self.k()}"

    def p(self) -> np.ndarray:
        """
        Access method for current position vector of the ray

        :return: np.ndarray, Ray position
        """
        return self._position[-1]

    def k(self) -> np.ndarray:
        """
        Access method for current direction vector of the ray

        :return: np.ndarray, Ray direction
        """
        return self._direction[-1]

    def append(self, p: np.ndarray | None, k: np.ndarray | None) -> None:
        """
        For a given p and k, appends these to hidden attributes _position and _direction. Either of these may be None
        type to terminate the ray

        :param p: p vector to append, or None to terminate
        :param k: k vector to append, or None to terminate
        :return: None
        """
        if p is None:
            self._position = np.append(self._position, [np.array([None, None, None])], axis=0)
        else:
            self._position = np.append(self._position, [p], axis=0)
        if k is None:
            self._direction = np.append(self._direction, [np.array([None, None, None])], axis=0)
        else:
            self._direction = np.append(self._direction, [k], axis=0)

    def vertices(self) -> np.ndarray:
        """
        Access method to the whole position vector, returning all points

        :return: np.ndarray, All points along the way
        """
        return self._position

    def is_terminated(self) -> bool:
        """
        Method to determine whether ray is terminated.

        A ray is terminated if any of p or k is None.

        If both p and K have None values in them, the ray is fully terminated after having attempted to intersect
        but point of intersection was not found. If p is defined (not None) but k is None, TIR has occurred so only the
        point of contact with the refracting object is given.

        :return: True if terminated, False otherwise
        """
        if None in self.p() or None in self.k():
            return True
        else:
            return False


class OpticalElement:
    def __init__(self):
        self._n1 = self._n2 = 1

    def normal(self, p_vector: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Must implement normal method")

    def _intercept(self, ray: Ray) -> None | np.ndarray:
        """
        Works out the nearest point of interception between the ray's origin and the Optical Element's surface. Returns
        array if this exists, None otherwise. Takes care of case where curvature = 0.

        :param Ray ray: Ray object to calculate the intersection of
        :return: np.ndarray object of the point of intersection in 3D space, None if this does not exist
        """
        def apt_check(intersection_pt):
            if vector_magnitude(intersection_pt - np.array([0, 0, intersection_pt[2]])) <= self._aperture:
                return intersection_pt
            else:
                return None

        p = ray.p()
        k = ray.k()
        k_hat = normalise(k)
        r = p - np.array([0, 0, self._z0])

        if self._curvature == 0:
            length = (- np.dot(r, z_hat)) / np.dot(k_hat, z_hat)
            return apt_check(p + length * k_hat)  # Checks that intersection is within aperture

        radius = 1 / self._curvature
        centre_lens = self._z0 + radius
        r = r + np.array([0, 0, radius])  # This is vector r as marked in the diagram before Task 4
        r = p - np.array([0, 0, centre_lens])
        r_dot_k = np.dot(r, k_hat)

        var_discriminant = r_dot_k ** 2 - (vector_magnitude(r) ** 2 - radius ** 2)

        if abs(var_discriminant) < 1e-7:
            # If there is only one intersection point i.e. quadratic nature is 0, Takes care of floating point operation
            # error
            return apt_check(p + r_dot_k * k_hat)  # in this case length = r_dot_k

        elif abs(var_discriminant) > 1e-7:  # There are 2 intersection points, discriminant is positive, so we can
            # safely sqrt. We select the smallest of the absolute value of the 2 lengths

            length_lst = [- r_dot_k + np.sqrt(var_discriminant), - r_dot_k - np.sqrt(var_discriminant)]

            if self._curvature > 0:
                length = min(length_lst)
            else:
                if abs(self._z0 - p[2]) < abs(radius):
                    length = length_lst[0]
                else:
                    length = max(length_lst)
            return apt_check(p + length * k_hat)

        #  Otherwise there is no intersection, so returns None

    def propagate_ray(self, ray: Ray) -> tuple:
        """
        Generalised propagate ray function to be properly implmented in sub-classes

        Raises an error if no intersection is found, or total internal reflection occurs, especially useful when
        debugging. This error does not break any functionality.

        :param Ray ray: ray to propagate
        :return: Tuple containing p and k vectors to append
        """

        new_p = self._intercept(ray)
        if new_p is None:  # Warns about the ray not intersecting, continues
            warnings.warn(f"\nNo intersection found for ray with p={ray.p()}, k={ray.k()}\nRay terminated")
            return None, None  # Return statement to exit

        normal_vec = self.normal(new_p)
        new_k = snell_refraction(normalise(ray.k()), normalise(normal_vec), self._n1, self._n2)
        if new_k is None:
            warnings.warn(f"\nTotal Internal Reflection for ray with p={ray.p()}, k={ray.k()}\nRay terminated")
            return new_p, None  # Return statement to exit

        return new_p, new_k


class SphericalRefraction(OpticalElement):

    def __init__(self, z0: float, curv: float, n1: float, n2: float, aperture: float):
        """
        Constructor for SphericalRefraction Class

        :param float z0: The intercept of the surface with the z-axis
        :param float curv: The curvature of the surface (1 / radius)
        :param float n1: Refractive index of non-spherical space
        :param float n2: Refractive index of sphere
        :param float aperture: The maximum extent of the surface from the optical axis
        :return: None
        """
        super().__init__()
        self._z0 = z0
        self._curvature = curv
        self._n1 = n1
        self._n2 = n2
        self._aperture = aperture

    def normal(self, p_vector: np.ndarray) -> np.ndarray:
        """
        Works out the normal of the surface at a given point (p_vector)

        :param np.ndarray p_vector: The point to work out the normal to the surface at
        :return: np.ndarray representing the normal
        """
        if self._curvature == 0:
            return np.array([0, 0, 1])
        return normalise(p_vector - (np.array([0, 0, self._z0 + 1 / self._curvature])))

    def propagate_ray(self, ray: Ray) -> None:
        """
        Inherits parent class' propagate_ray method propagating the ray after having intersected with the spherical
        refractive object in question (this object).

        Calculates the intersection point (if exists), appends this point to the ray, and using global snell_refraction
        function works out the direction of the ray and appends this too.

        Terminates ray if there is no intersection point or if TIR occurs. See Ray.is_terminated() docs for further info
        on ray termination.

        :param Ray ray: Ray object to propagate
        :return: None
        """
        if ray.is_terminated():
            warnings.warn(f"\nRay {ray} is terminated, cannot propagate")
            return None
        pos, direct = super().propagate_ray(ray)
        ray.append(pos, direct)


class OutputPlane(OpticalElement):

    def __init__(self, z0: float, aperture: float) -> None:
        super().__init__()
        self._z0 = z0
        self._aperture = aperture
        self._curvature = 0

    def normal(self, p_vector: np.ndarray) -> np.ndarray:
        """
        Works out the normal of the surface at a given point (p_vector)

        :param np.ndarray p_vector: The point to work out the normal to the surface at
        :return: np.ndarray representing the normal
        """
        return np.array([0, 0, -1])

    def propagate_ray(self, ray: Ray) -> None:
        """
        Inherits parent class' propagate_ray method propagating the ray after having intersected with the output plane
        in question (this object).

        Calculates the intersection point (if exists), appends this point to the ray, and terminates the ray by
        appending None to the ray's direction vector. See Ray.is_terminated() docs for further info on ray termination.

        Its calling will have no effect if the ray is already terminated.

        :param Ray ray: Ray object to propagate
        :return: None
        """
        if ray.is_terminated():
            warnings.warn(f"\nRay {ray} is terminated, cannot propagate")
            return None
        pos, _ = super().propagate_ray(ray)
        ray.append(pos, None)  # Marks ray as terminated


class Lens:
    def __init__(self, z0: float, curv_in: float, curv_out: float, n_lens, aperture: float, n_out: float = 1):
        """
        Class constructor.

        Works out the point of intersection of the right-most refractive plane automatically.

        :param z0: Intersection of the initial refractive plane with the z-axis
        :param curv_in: Curvature of initial refractive plane
        :param curv_out: Curvature of final refractive plane
        :param n_lens: Refractive index of lens
        :param aperture: Aperture, or size of the lens (radius)
        :param n_out: Refractive index outside the lens (air or vacuum generally)
        :raise ValueError: If both sides of the lens have 0 curvature
        """
        if curv_in == curv_out == 0:
            raise ValueError("Both sides of the lens cannot be of 0 curvature")

        self._z0_in = z0
        self._curv_in = curv_in
        self._n_lens = n_lens
        self._n_air = n_out
        self._aperture = aperture
        self._curv_out = curv_out
        self._z0_out = None  # This is defined in the method z0_out_generator

        self.z0_out_generator()

        self.surface_in = SphericalRefraction(self._z0_in, self._curv_in, self._n_air, self._n_lens, self._aperture)

        self.surface_out = SphericalRefraction(self._z0_out, self._curv_out, self._n_air, self._n_lens, self._aperture)

    def z0_out_generator(self):
        """
        Works out z-axis intersection of the final refractive plane

        :return: None
        :raise NotImplementedError: if not implemented by sub-class
        """
        raise NotImplementedError

    def surface_thickness(self, curvature: float) -> float:
        """
        Useful tool to work out the thickness of a refractive surface given its curvature and the lens' aperture
        (radius)

        :param float curvature: Curvature of the refractive surface to work out the thickness of
        :return: float - Thickness of refractive surface
        """
        return (1 / abs(curvature)) - \
               (1 / abs(curvature)) * np.sin(np.arccos(self._aperture * abs(curvature)))

    def _propagate(self, ray: Ray) -> None:
        """
        Propagates ray through the lens in question

        :param Ray ray: Ray to propagate
        :return: None
        """
        self.surface_in.propagate_ray(ray)
        self.surface_out.propagate_ray(ray)

    def plot_xz(self, ray_lst: list, outpt_pln: OutputPlane | None = None) -> None:
        """
        Plots graph in the xz plane (ray traces).

        Must use plt.show() statement after the use of this method, this is useful to further customise the plot (adding
        title etc.)

        :param list ray_lst: List containing ray objects to plot
        :param outpt_pln: Optional, Output plane if desired
        :return: None
        """
        for ray in ray_lst:
            self._propagate(ray)
            if outpt_pln is not None:
                outpt_pln.propagate_ray(ray)
            points = ray.vertices()
            plt.plot([i[2] for i in points],
                     [i[0] for i in points], '-x', color='k') if None not in points else None

    def plot_xy(self, ray_lst: list, outpt_pln: OutputPlane) -> None:
        """
        Plots scatter graph in the z plane (dots)

        Must use plt.show() statement after the use of this method, this is useful to further customise the plot (adding
        title etc.)

        :param ray_lst: List containing ray objects to plot
        :param outpt_pln: Output plane to plot the graph on. Unlike in plot_xz, this is mandatory here
        :return: None
        """
        plt.figure(figsize=(6, 6))
        for ray in ray_lst:
            self._propagate(ray)
            outpt_pln.propagate_ray(ray)
            points = ray.p()
            plt.scatter(points[0], points[1], color='k')


class PlanoConvex(Lens):
    def __init__(self, z0: float, curv_in: float, curv_out: float, n_lens:float, aperture: float,n_out: float = 1):
        """
        Class constructor.

        Works out the point of intersection of the right-most refractive plane automatically.

        :param z0: Intersection of the initial refractive plane with the z-axis
        :param curv_in: Curvature of initial refractive plane
        :param curv_out: Curvature of final refractive plane
        :param n_lens: Refractive index of lens
        :param aperture: Aperture, or size of the lens (radius)
        :param n_out: Refractive index outside the lens (air or vacuum generally)
        :raise ValueError: If input parameters for curvature do not make a plano-convex lens
        """
        if curv_in == 0 > curv_out or curv_in > 0 == curv_out:
            Lens.__init__(self, z0, curv_in, curv_out, n_lens, aperture, n_out)
        else:
            raise ValueError("The input parameters for curvature do not make a Plano-Convex Lens")

    def z0_out_generator(self):
        """
        Works out z-axis intersection of the final refractive plane for a plano-convex configuration.

        :return: None
        """
        if self._curv_in == 0:
            self._z0_out = self._z0_in + self.surface_thickness(self._curv_out)

        else:
            self._z0_out = self._z0_in + self.surface_thickness(self._curv_in)


class BiConvex(Lens):
    def __init__(self, z0: float, curv_in: float, curv_out: float, n_lens: float, aperture: float, n_out: float = 1):
        """
        Class constructor.

        Works out the point of intersection of the right-most refractive plane automatically.

        :param z0: Intersection of the initial refractive plane with the z-axis
        :param curv_in: Curvature of initial refractive plane
        :param curv_out: Curvature of final refractive plane
        :param n_lens: Refractive index of lens
        :param aperture: Aperture, or size of the lens (radius)
        :param n_out: Refractive index outside the lens (air or vacuum generally)
        :raise ValueError: If input parameters for curvature do not make a biconvex lens
        """
        if curv_in > 0 > curv_out:
            Lens.__init__(self, z0, curv_in, curv_out, n_lens, aperture, n_out)
        else:
            raise ValueError("The input parameters for curvature do not make a Biconvex Lens")

    def z0_out_generator(self):
        """
        Works out z-axis intersection of the final refractive plane for a biconvex configuration.

        :return: None
        """
        self._z0_out = self._z0_in + 2 * self.surface_thickness(self._curv_in)


def snell_refraction(inc_v: np.ndarray, surf_n: np.ndarray, n_1: float, n_2: float):
    """
    Computational method of for Snell's law. With given incident ray direction vector, refractive object normal vector,
    incident refractive index, refracting material refractive index it will output a unit directional vector of the
    refracted ray. Returns none if Total Internal Reflection occurs

    :param np.ndarray inc_v: Incident ray direction vector
    :param np.ndarray surf_n: Plane normal vector
    :param n_1: Refractive index of incident material
    :param n_2: Refractive index of refracting material
    :return: Direction vector of refracted ray, None if TIR occurs
    """
    if 0 in (n_1, n_2):
        raise ValueError("Refractive indexes cannot be zero")
    surf_n_nor = normalise(surf_n)
    inc_v_nor = normalise(inc_v)
    theta1 = np.arcsin(vector_magnitude(np.cross(inc_v_nor, surf_n_nor)))
    if np.sin(theta1) > (n_2 / n_1):  # Checks for TIR
        return None

    theta2 = np.arcsin(n_1 / n_2 * vector_magnitude(np.cross(inc_v_nor, surf_n_nor)))
    ref_v = (n_1 / n_2) * inc_v_nor + ((n_1 / n_2) * np.cos(theta1) - np.cos(theta2)) * surf_n_nor
    return normalise(ref_v)


def vector_magnitude(vector: np.ndarray) -> float:
    """
    Works out the vector scalar magnitude

    :param np.ndarray vector: The vector to find the magnitude of
    :return: Vector magnitude
    """
    return np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)


def normalise(dir_vect: np.ndarray) -> np.ndarray:
    """
    Returns normalised direction vector

    :param np.ndarray dir_vect: Input direction vector to normalise
    :return: Normalised direction vector
    """
    return dir_vect / vector_magnitude(dir_vect)


def col_generator(z0: float = 0, rad: float = 5, dist_pts: float = 0.5, k_vec: None | np.ndarray = None) -> list:
    """
    Generates a collimated bundle of equally spaced rays travelling in the z-direction by default using polar
    coordinates.

    Although it is theoretically impossible to make an equally spaced bundle of rays, by increasing the number of rays
    by 6 ever new ring we can get the closest approximation

    Calling the function without parsing any parameters will output a default ray bundle of radius 5

    :param float z0: the z point to start the bundle of rays at
    :param float rad: the radius of the ray bundle
    :param float dist_pts: distance between points
    :param None | np.ndarray k_vec: direction vector of the rays (defaults to z_hat)
    :return: list - containing all individual ray objects
    """
    if k_vec is None:  # Properly dealing with mutable objects passed as kwargs
        k = np.array([0, 0, 1])
    else:
        k = k_vec
    ray_lst = [Ray(np.array([0, 0, z0]), k)]

    pts_incr = 6  # 6 is the best approximation to "uniformly distributed". This is impossible in reality
    num_of_shells = int(rad / dist_pts)
    for i in range(1, num_of_shells):
        num_pts = i * pts_incr  # The number of rays in the "shell"
        radius = i * dist_pts
        ray_lst += (Ray(pos=np.array([radius * np.cos(2 * np.pi * dot / num_pts),
                                      radius * np.sin(2 * np.pi * dot / num_pts),
                                      z0]),
                        direc=k)
                    for dot in range(num_pts))

    return ray_lst


def rms(ray_lst: list[Ray]) -> float:
    """
    Works out the RMS radius of a list of ray objects

    :param list[Ray] ray_lst: List of rays to calculate the RMS radius of
    :return: float - RMS radius of given rays
    """
    squares_list = []
    for ray in ray_lst:
        if None not in ray.p():
            squares_list.append(ray.p()[0] ** 2 + ray.p()[1] ** 2)

    return np.sqrt(np.mean(squares_list))


def rms_plotter(ray_lst: list[Ray], arange: np.ndarray) -> tuple:
    """
    Helper function to easily debug and further expand program. Plots rms values of rays over a range of previously
    given values.

    Helps find the point of minimum RMS (focal point)

    :param ray_lst: List of rays to calculate the RMS radius of
    :param np.ndarray arange: range of z-values to iterate over
    :return: Tuple containing z-values and their corresponding RMS values
    """
    lst_of_rays = copy.deepcopy(ray_lst)
    z_vals = []
    rms_vals = []
    for z_val in arange:
        ray_set = copy.deepcopy(lst_of_rays)
        out_pln = OutputPlane(z_val, 50)
        for ray in ray_set:
            out_pln.propagate_ray(ray)
        rms_vals.append(rms(ray_set))
        z_vals.append(z_val)

    return z_vals, rms_vals
