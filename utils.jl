using LinearAlgebra, BlockDiagonals

# Canonical bases in R^3
e(i) = I[1:3, i]

# Trigonometric vector 
u(angle) = [sin(angle); cos(angle)]
angle(u) = atan(u[1] / u[2])

# Rotation matrix in R^2
R0(u) = [
    u[2] -u[1]
    u[1] u[2]
]

# Rotation matrix in R^3, around West-East axis
R1(u) = BlockDiagonal([ones(1, 1), R0(u)])

# Rotation matrix in R^3, around Up axis
R3(u) = BlockDiagonal([R0(u), ones(1, 1)])

# Panel vertices position when lying flat on the ground, 
# with the center corresponding to the origin
base_vertices(width, depth) = [
    -width width width -width
    -depth -depth depth depth
    0 0 0 0
]

vertices(width, depth, center, planar_azimuth_vec, tilt_vec) =
    R3(planar_azimuth_vec)' * R1(tilt_vec) * base_vertices(width, depth) .+ center

light_beam(solar_azimuth_vec, solar_elevation_vec) =
    -R3(solar_azimuth_vec)' * R1(solar_elevation_vec) * e(2)

discrete_sin_wave(n) = sin(π * n / 2)

# Irradiance terms
proj_incidence_vec(panel, sun) =
    [cos(sun.elevation) * cos(sun.azimuth - panel.azimuth), sin(sun.elevation)]

irr_normal(panel, sun, tilt_vec) = sun.dni * proj_incidence_vec(panel, sun)' * tilt_vec
irr_diff(panel, sun, tilt_vec) = 0.5sun.dni * (1 + tilt_vec[2])
irr_ground(panel, sun, albedo, tilt_vec) = 0.5albedo * sun.ghi * (1 - tilt_vec[2])

area(panel) = panel.width * panel.depth

panel_power(panel, sun, albedo, tilt_vec) =
    area(panel) * (
        irr_normal(panel, sun, tilt_vec) +
        irr_diff(panel, sun, tilt_vec) +
        irr_ground(panel, sun, albedo, tilt_vec)
    )

# Shadow convex hull linear components
function a(panel, crop, sun, vertex_index)
    a1 =
        0.5panel.width *
        panel.depth *
        [cos(sun.azimuth - panel.azimuth) / tan(sun.elevation), 1]
    a21 =
        cos(sun.azimuth) / tan(sun.elevation) * (panel.pos_east - crop.pos_east) +
        sin(sun.azimuth) / tan(sun.elevation) * (crop.pos_north - panel.pos_north)
    a22 =
        cos(panel.azimuth) * (panel.pos_east - crop.pos_east) +
        sin(panel.azimuth) * (crop.pos_north + panel.pos_north) -
        panel.height * sin(sun.azimuth - panel.azimuth) / tan(sun.elevation)
    return a1 + discrete_sin_wave(vertex_index - 1) * panel.depth * [a21, a22]
end

function b(panel, crop, sun, vertex_index)
    b1 = panel.height * cos(sun.azimuth - panel.azimuth) / tan(sun.elevation)
    b2 =
        sin(panel.azimuth) * (crop.pos_east - panel.pos_east) +
        cos(panel.azimuth) * (crop.pos_north - panel.pos_north)
    return panel.width * discrete_sin_wave(vertex_index) * (b1 + b2)
end

# If the returned expression is nonnegative, the crop is inside the 
# halfplane identified by the vertex_index
shadowing_condition(panel, crop, sun, vertex_index, tilt_vec) =
    a(panel, crop, sun, vertex_index)' * tilt_vec + b(panel, crop, sun, vertex_index)

# If the returned expression is nonnegative, the sun is hitting the panel
# on the side there the PV cells are mounte, i.e., the direct irradiance 
# flux is nonnegative.
proj_incidence_condition(panel, sun, tilt_vec) = proj_incidence_vec(panel, sun)' * tilt_vec

# Check (probabilistically) whether a panel can shadow a crop, 
# given the sun position
can_panel_shadow_crop(panel, crop, sun, num_samples) = any(
    angle -> (
        all(j -> shadowing_condition(panel, crop, sun, j, u(angle)) >= 0, 1:4) &
        (proj_incidence_condition(panel, sun, u(angle)) >= 0)
    ),
    range(start = -π / 2, stop = π / 2, length = num_samples),
)

# big-M and small-M for logical to integer constraint
# reformulation
function big_M(panel, crop, sun, vertex_index) 
  _a = a(panel, crop, sun, vertex_index)
  return - b(panel, crop, sun, vertex_index) + _a' * sign.(_a)
end

function small_M(panel, crop, sun, vertex_index) 
  _a = a(panel, crop, sun, vertex_index)
  return - b(panel, crop, sun, vertex_index) - _a' * sign.(_a)
end
