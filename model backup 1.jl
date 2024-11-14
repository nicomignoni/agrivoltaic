using LinearAlgebra, BlockDiagonals, JuMP, MosekTools

struct Panel
    width::Real
    depth::Real
    pos::Vector
    azimuth::Real
end

struct Crop
    pos::Vector
    to_shadow::Bool
end

struct Sun
    time::String
    dni::Real
    dhi::Real
    ghi::Real
    azimuth::Real
    elevation::Real
end

# --------------------- Panel ---------------------

panel_area(panel::Panel) = panel.width * panel.depth

panel_power(panel::Panel, sun::Sun, albedo::Real, tilt_vec) =
    panel_area(panel) * (
        irr_normal(panel, sun, tilt_vec) +
        irr_diff(sun, tilt_vec) +
        irr_ground(sun, albedo, tilt_vec)
    )

base_vertices(panel::Panel) = [
    -panel.width panel.width panel.width -panel.width
    -panel.depth -panel.depth panel.depth panel.depth
    0 0 0 0
]

vertices(panel::Panel, tilt_vec) =
    R3(u(panel.azimuth))' * R1(tilt_vec) * base_vertices(panel) .+ panel.pos

# --------------------- Sun & Irradiance ---------------------

light_beam(sun::Sun) = -R3(u(sun.azimuth))' * R1(u(sun.elevation)) * e(2)

proj_incidence_vec(panel, sun) =
    [cos(sun.elevation) * cos(sun.azimuth - panel.azimuth), sin(sun.elevation)]

irr_normal(panel::Panel, sun::Sun, tilt_vec) = sun.dni * proj_incidence_vec(panel, sun)' * tilt_vec

irr_diff(sun::Sun, tilt_vec) = 0.5sun.dni * (1 + tilt_vec[2])

irr_ground(sun::Sun, albedo::Real, tilt_vec) = 0.5albedo * sun.ghi * (1 - tilt_vec[2])

# --------------------- Geometry --------------------

raw"Canonical bases in $\mathbf{R}^3$."
e(i::Int) = I[1:3, i]

raw"Trigonometric vector."
u(angle::Real) = [sin(angle); cos(angle)]

raw"Angle represented by the trigonometric vector"
angle(u::Vector) = atan(u[1] / u[2])

raw"Rotation matrix in $\mathbf{R}^2.$"
R0(u::Vector) = [
    u[2] -u[1]
    u[1] u[2]
]

raw"Rotation matrix in $\mathbf{R}^3$, around West-East (i.e., X) axis."
R1(u::Vector) = BlockDiagonal([ones(1, 1), R0(u)])

raw"Rotation matrix in $\mathbf{R}^3$, around Up axis."
R3(u::Vector) = BlockDiagonal([R0(u), ones(1, 1)])

raw"""
Panel vertices position when lying flat on the ground, 
with the center corresponding to the origin.
"""
discrete_sin_wave(n::Int) = sin(π * n / 2)

# --------------------- Control problem ---------------------

# Shadow convex hull linear components
function a(panel::Panal, crop::Crop, sun::Sun, vertex_index)
    a1 =
        0.5panel.width *
        panel.depth *
        [cos(sun.azimuth - panel.azimuth) / tan(sun.elevation), 1]
    a21 =
        cos(sun.azimuth) / tan(sun.elevation) * (panel.pos[1] - crop.pos[1]) +
        sin(sun.azimuth) / tan(sun.elevation) * (crop.pos[2] - panel.pos[2])
    a22 =
        cos(panel.azimuth) * (panel.pos[1] - crop.pos[1]) +
        sin(panel.azimuth) * (crop.pos[2] + panel.pos[2]) -
        panel.pos[3] * sin(sun.azimuth - panel.azimuth) / tan(sun.elevation)
    return a1 + discrete_sin_wave(vertex_index - 1) * panel.depth * [a21, a22]
end

function b(panel, crop, sun, vertex_index)
    b1 = panel.pos[3] * cos(sun.azimuth - panel.azimuth) / tan(sun.elevation)
    b2 =
        sin(panel.azimuth) * (crop.pos[1] - panel.pos[1]) +
        cos(panel.azimuth) * (crop.pos[2] - panel.pos[2])
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
    return -b(panel, crop, sun, vertex_index) + _a' * sign.(_a)
end

function small_M(panel, crop, sun, vertex_index)
    _a = a(panel, crop, sun, vertex_index)
    return -b(panel, crop, sun, vertex_index) - _a' * sign.(_a)
end

# Optimal control open-loop
function control(sun, panels, crops, albedo, gamma, num_samples = 10, epsilon = 1e-4)
    print("[Sun @ $(sun.time)] Start constructing problem... ")
    prob = Model(MosekTools.Optimizer)
    set_attribute(prob, "QUIET", true)

    # Consider a pair (panel, crop) only if the panel can actually
    # shadow the crop
    panel_crop_pairs = [
        (i, k) for (i, panel) in enumerate(panels) for (k, crop) in enumerate(crops) if
        can_panel_shadow_crop(panel, crop, sun, num_samples)
    ]

    # Planar tilt vector, with bounds
    @variable(prob, z - 2 <= tilt_vec[1:length(panels), z = 1:2] <= 1)

    # Is the crop position shadowed?
    @variable(prob, is_crop_shadowed[1:length(crops)], Bin, start = 0)

    # Fundamental trigonometric equality constraint
    if gamma > 0
        # Relaxed form
        @constraint(prob, tilt_vec[:, 1] .^ 2 .+ tilt_vec[:, 2] .^ 2 .<= 1)
    end

    # Avoid light to hit the back of the panel 
    # (negative irradiance flux -> negative power) 
    @constraint(
        prob,
        [(i, panel) in enumerate(panels)],
        proj_incidence_condition(panel, sun, tilt_vec[i, :]) >= 0
    )

    num_modules_shadowing_crop = zeros(AffExpr, length(crops))

    # Does the panel shadow the crop position?
    does_module_shadow_crop = Dict{NTuple{2,Int},VariableRef}()

    # Does the panel's halfplane shadows the crop position?
    does_module_hp_shadow_crop = Dict{NTuple{3,Int},VariableRef}()

    for (i, k) in panel_crop_pairs
        does_module_shadow_crop[(i, k)] = @variable(prob, binary = true)
        num_modules_shadowing_crop[k] += does_module_shadow_crop[(i, k)]

        panel, crop = panels[i], crops[k]

        # Shadow indicators definitions' constraints
        for j = 1:4
            does_module_hp_shadow_crop[(i, j, k)] = @variable(prob, binary = true)
            @constraints(
                prob,
                begin
                    -shadowing_condition(panel, crop, sun, j, tilt_vec[i, :]) <=
                    big_M(panel, crop, sun, j) *
                    (1 - does_module_hp_shadow_crop[(i, j, k)])
                    -shadowing_condition(panel, crop, sun, j, tilt_vec[i, :]) >=
                    epsilon +
                    (small_M(panel, crop, sun, j) - epsilon) *
                    does_module_hp_shadow_crop[(i, j, k)]
                end
            )
        end

        # Is the crop contained in all the modules' halfplanes (i.e., is 
        # the module shadowing the crop)? 
        @constraints(
            prob,
            begin
                sum(j -> does_module_hp_shadow_crop[(i, j, k)], 1:4) <=
                does_module_shadow_crop[(i, k)] + 3
                sum(j -> does_module_hp_shadow_crop[(i, j, k)], 1:4) >=
                4does_module_shadow_crop[(i, k)]
            end
        )
    end

    # Does any module shadow the crop?
    @constraints(prob, begin
        num_modules_shadowing_crop .<= length(panels) * is_crop_shadowed
        num_modules_shadowing_crop .>= is_crop_shadowed
    end)

    # Objective function terms
    total_power = sum(
        panel_power(panel, sun, albedo, tilt_vec[i, :]) for (i, panel) in enumerate(panels)
    )

    # In JuMP, the l1 norm is explicitly defined throught the canonical cone formulation
    shadow_reference = stack(crop.to_shadow for crop in crops)
    @variable(prob, reference_distance)
    @constraint(
        prob,
        [reference_distance; is_crop_shadowed - shadow_reference] in
        MOI.NormOneCone(1 + length(crops))
    )
    print("Done. ")

    # Solve
    @objective(prob, Max, gamma * total_power - (1 - gamma) * reference_distance)
    optimize!(prob)
    println("Status: $(termination_status(prob)).")

    return value.(u), value.(is_crop_shadowed)
end
