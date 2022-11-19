'''
:Date: 17-11-2022
:Author: Daniel Weronski Falco, Blackett Laboratory

First module

Classes:


Functions:

'''

import numpy as np
import warnings
from typing import Union


class Ray:
    def __init__(self, pos: np.ndarray, direc: np.ndarray) -> None:
        self._position = np.array([pos])
        self._direction = np.array([direc])

    def p(self) -> np.ndarray:
        """
        Access method for current position vector
        :return: np.ndarray, Ray position
        """
        return self._position[-1]

    def k(self) -> np.ndarray:
        """
        Access method for current direction vector
        :return: np.ndarray, Ray direction
        """
        return self._direction[-1]

    def append(self, p: Union[np.ndarray, None], k: Union[np.ndarray, None]) -> None:
        """
        For a given p and k, appends these to hidden attributes _position and _direction. Either of these may be None
        type to terminate the ray

        :param p: p vector to append, or None to terminate
        :param k: k vector to append, or None to terminate
        :return: None
        """
        self._position = np.append(self._position, [p], axis=0)
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
        A ray is terminated if any of p or k is None. It may be that p is not None but k is, meaning TIR has occurred

        :return: True if terminated, False otherwise
        """
        if None in self.p() or None in self.k():
            return True
        else:
            return False


class OpticalElement:
    def __init__(self):
        pass

    # def propagate_ray(self, ray):
    #     "propagate a ray through the optical element"


class SphericalRefraction:  # (OpticalElement):

    def __init__(self, z0: float, curv: float, n1: float, n2: float, a_radius: float) -> None:
        """
        --

        :param float z0: The intercept of the surface with the z-axis
        :param float curv: The radius (of curvature) of the surface
        :param float n1: Refractive index of non-spherical space
        :param float n2: Refractive index of sphere
        :param float a_radius: The maximum extent of the surface from the optical axis
        :return: None
        """
        self._z0 = z0
        self._curvature = curv
        self._n1 = n1
        self._n2 = n2
        self.aperture_radius = a_radius

    def intercept(self, ray: Ray):
        """
        Works out the nearest point of interception between the ray's origin and the Optical Element's surface. Returns
        array if this exists, None otherwise.

        :param Ray ray: Ray object to calculate the intersection of
        :return: np.ndarray object of the point of intersection in 3D space, None if this does not exist
        """

        def discriminant() -> float:
            """
            :return: Value of b^2 - 4ac (a=1), in the context of the quadratic equation to find out the length between
            p and the point of intersection
            """
            return np.dot(r, k_hat) ** 2 - (vector_magnitude(r) ** 2 - self._curvature ** 2)

        p = ray.p()
        k = ray.k()
        k_hat = normalise(k)
        r = np.array([0, 0, self._z0]) - p
        r_dot_k = np.dot(r, k_hat)

        if discriminant() < -1e-7:  # There is no intersection, takes care of floating point operation error
            return None

        elif abs(discriminant()) < 1e-7:
            # If there is only one intersection point i.e. quadratic nature is 0,Takes care of floating point operation
            # error

            return p + r_dot_k * k_hat  # in this case length = r_dot_k

        else:  # There are 2 intersection points
            # There can be 2 lengths of which we must select the smallest absolute value
            length_lst = [abs(- r_dot_k + np.sqrt(discriminant())), abs(- r_dot_k - np.sqrt(discriminant()))]
            length = min(length_lst)
            return p + length * k_hat

    def propagate_ray(self, ray: Ray) -> None:
        """
        # TODO: Make sure to properly document the None returns
        # TODO: Make sure to deal with the ray termination stuff


        :param Ray ray: Ray object to propagate
        :return: None
        """
        new_p = self.intercept(ray)
        if new_p is None:  # Warns about the ray not intersecting, continues
            ray.append(np.array([None, None, None]), np.array([None, None, None]))
            warnings.warn(f"No intersection found for ray with p={ray.p()}, k={ray.k()}\nRay terminated")
            return None  # Return statement to exit

        normal_vec = np.array([new_p[0], new_p[1], new_p[2] - self._z0])
        new_k = snell_refraction(normalise(ray.k()), normalise(normal_vec), self._n1, self._n2)
        if new_k is None:
            ray.append(new_p, np.array([None, None, None]))
            warnings.warn(f"Total Internal Reflection for ray with p={ray.p()}, k={ray.k()}\nRay terminated")
            return None  # Return statement to exit
        
        ray.append(new_p, new_k)


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
    if not all([n_1, n_2]):
        raise ValueError("Refractive indexes cannot be zero")

    theta1 = np.arccos(np.dot(inc_v, surf_n) / (vector_magnitude(inc_v) * vector_magnitude(surf_n)))
    if np.sin(theta1) > (n_2 / n_1):  # Checks for TIR
        return None

    theta2 = np.arcsin(n_1 * np.sin(theta1) / n_2)
    ref_v = (n_1 / n_2) * inc_v + ((n_1 / n_2) * np.cos(theta1) - np.cos(theta2)) * surf_n
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
    return np.array([dir_vect[0] / vector_magnitude(dir_vect),
                     dir_vect[1] / vector_magnitude(dir_vect),
                     dir_vect[2] / vector_magnitude(dir_vect)])



########################################################################################################################
# snell_refraction(np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0]), np.array([1, 0, 0]), 1.0, 1.5)

ray1 = Ray([10, 0, 0], [-2, 0, 0])

# ray1.append(np.array([1, 2, 3]), np.array([1, 2, 3]))

sphere = SphericalRefraction(0, 5, 1, 1.5, 0)
print(sphere.intercept(ray1))
sphere.propagate_ray(ray1)
print(ray1.is_terminated())
# print()
# ########################################################################################################################
#
# print(sphere.intercept(ray1))
# print(type(thing.get_position()))
