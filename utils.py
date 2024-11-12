import numpy as np
import scipy as sp

# Canonical base
e = lambda i: np.eye(1,3,i).T

# Trigonometric vector
u = lambda angle: np.array([np.sin(angle), np.cos(angle)])

# Rotation matrix in R^2
R0 = lambda u: np.block([
    [u[1], -u[0]],
    [u[0],  u[1]]
])

# Rotation matrix in R^3, around East-West axis
R1 = lambda u: sp.linalg.block_diag(1, R0(u))

# Rotation matrix in R^3, around Up axis
R3 = lambda u: sp.linalg.block_diag(R0(u), 1)

H0 = R0(u(np.pi/2))
H3 = R3(u(np.pi/2))

def base_vertices(width, depth): 
    return 0.5*np.block([
        [-width,  width, width, -width],
        [-depth, -depth, depth,  depth],
        [     0,      0,     0,      0]
    ])

def vertices(width, depth, center, planar_azimuth_vec, planar_tilt_vec):
    return R3(planar_azimuth_vec).T @ R1(planar_tilt_vec).T @ base_vertices(width, depth) + center

def light_beam(solar_azimuth_vec, solar_elevation_vec):
    return -R3(solar_azimuth_vec).T @ R1(solar_elevation_vec) @ e(1)

def projection_matrix(solar_azimuth_vec, solar_elevation_vec):
    beam = light_beam(solar_azimuth_vec, solar_elevation_vec)
    return np.eye(3) - beam @ e(2).T / beam[2]

# # Normal irradiance component
# def projected_incidence(planar_azimuth: float, planar_tilt: float, solar_azimuth: float, solar_elevation: float):
#     return np.sin(planar_tilt)*np.cos(solar_elevation)*np.cos(solar_azimuth - planar_azimuth) + np.sin(solar_elevation)*np.cos(planar_tilt)
#
# @multimethod
# def projected_incidence(planar_azimuth_u: ExtendedIterable, planar_tilt_u: ExtendedIterable, solar_azimuth: ExtendedIterable, solar_elevation_u: ExtendedIterable):
#     return planar_tilt_u[0]*solar_elevation_u[1]*(solar_azimuth_u[0]*planar_azimuth_u[0] + solar_azimuth_u[1]*planar_azimuth_u[1]) + solar_elevation_u[0]*planar_tilt_u[1]
#
# def irradiance_normal(dni, planar_azimuth, planar_tilt, solar_azimuth, solar_elevation):
#     return dni*projected_incidence(planar_azimuth, planar_tilt, solar_azimuth, solar_elevation)
#
# # Diffuse irradiance component
# @multimethod
# def irradiance_diffuse(dhi, planar_azimuth: float):
#     return 0.5*dni*(1 + np.cos(planar_azimuth))
#
# @multimethod
# def irradiance_diffuse(dhi, planar_azimuth_u: ExtendedIterable):
#     return 0.5*dhi*(1 + planar_azimuth_u[1])
#
# # Ground reflected irradiance component
# @multimethod
# def irradiance_ground(ghi, planar_azimuth: float):
#     return 0.5*ghi*(1 - np.cos(planar_azimuth))
#
# @multimethod
# def irradiance_ground(ghi, planar_azimuth_u: ExtendedIterable):
#     return 0.5*ghi*(1 - planar_azimuth[1])
#
