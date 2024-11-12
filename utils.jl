using LinearAlgebra, BlockDiagonals

# Canonical bases in R^3
e(i) = I[1:3, i]

# Trigonometric vector 
u(angle) = [sin(angle); cos(angle)]
angle(u) = atan(u[1] / u[2])

# Rotation matrix in R^2
R0(u) = [
  u[2] -u[1]
  u[1]  u[2]
]

# Rotation matrix in R^3, around West-East axis
R1(u) = BlockDiagonal([ones(1,1), R0(u)])

# Rotation matrix in R^3, around Up axis
R3(u) = BlockDiagonal([R0(u), ones(1,1)])

base_vertices(width, depth) = [
  -width  width  width -width
  -depth -depth  depth  depth
       0      0      0      0
]

vertices(width, depth, center, planar_azimuth_vec, planar_tilt_vec) = 
  R3(planar_azimuth_vec)' * R1(planar_tilt_vec) * base_vertices(width, depth) .+ center

light_beam(solar_azimuth_vec, solar_elevationation_vec) = 
  -R3(solar_azimuth_vec)' * R1(solar_elevationation_vec) * e(2)

discrete_sin_wave(n) = sin(π * n / 2)

# Irradiance
proj_incidence_vec(panel, sun) =
    [cos(sun.elevation) * cos(sun.azimuth - panel.azimuth), sin(sun.elevation)]

irr_normal(panel, sun, tilt_vec) = sun.dni * proj_incidence_vec(panel, sun)' * tilt_vec
irr_diff(panel, sun, tilt_vec) = 0.5sun.dni * (1 + tilt_vec[2])
irr_ground(panel, sun, albedo, tilt_vec) = 0.5albedo * sun.ghi * (1 - tilt_vec[2])

# Power
area(panel) = panel.width * panel.depth
panel_power(panel, sun, albedo, tilt_vec) =
    area(panel) * (
        irr_normal(panel, sun, tilt_vec) +
        irr_diff(panel, sun, tilt_vec) +
        irr_ground(panel, sun, albedo, tilt_vec)
    )
