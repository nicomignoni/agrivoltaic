using LinearAlgebra, BlockDiagonals, JuMP, SCIP

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

"""Panel's area"""
A(panel::Panel) = panel.width * panel.depth

"""Panel's absorbed power"""
p(panel::Panel, sun::Sun, albedo::Real, tilt_vec) =
    A(panel) * (
        irr_normal(panel, sun, tilt_vec) +
        irr_diff(sun, tilt_vec) +
        irr_ground(sun, albedo, tilt_vec)
    )

"""
Panel's vertices position when lying flat on the ground, 
with the center corresponding to the origin.
"""
v₀(panel::Panel) = [
    -panel.width panel.width panel.width -panel.width
    -panel.depth -panel.depth panel.depth panel.depth
    0 0 0 0
]

raw"""Panel's vertices position in $\mathbb{R}^3$$"""
v(panel::Panel, tilt_vec) = R₃(u(panel.azimuth))' * R₁(tilt_vec) * v₀(panel) .+ panel.pos

# --------------------- Sun & Irradiance ---------------------

# If the returned expression is nonnegative, the sun is hitting the panel
# on the side there the PV cells are mounted, i.e., the direct irradiance
# flux is nonnegative.
proj_incidence_condition(panel::Panel, sun::Sun, tilt_vec) =
    proj_incidence_vec(panel, sun)' * tilt_vec

light_beam(sun::Sun) = -R₃(u(sun.azimuth))' * R₁(u(sun.elevation)) * e(2)

proj_incidence_vec(panel::Panel, sun::Sun) =
    [cos(sun.elevation) * cos(sun.azimuth - panel.azimuth), sin(sun.elevation)]

irr_normal(panel::Panel, sun::Sun, tilt_vec) = sun.dni * proj_incidence_vec(panel, sun)' * tilt_vec

irr_diff(sun::Sun, tilt_vec) = 0.5sun.dni * (1 + tilt_vec[2])

irr_ground(sun::Sun, albedo::Real, tilt_vec) = 0.5albedo * sun.ghi * (1 - tilt_vec[2])

# --------------------- Geometry --------------------

raw"Canonical bases in $\mathbf{R}^3$."
e(i::Int) = I[1:3, i]

raw"Trigonometric vector."
u(θ::Real) = [sin(θ); cos(θ)]

raw"Angle represented by the trigonometric vector"
angle(u::Vector) = atan(u[1] / u[2])

raw"Rotation matrix in $\mathbf{R}^2.$"
R₀(u::Vector) = [
    u[2] -u[1]
    u[1] u[2]
]

raw"Rotation matrix in $\mathbf{R}^3$, around West-East (i.e., X) axis."
R₁(u::Vector) = BlockDiagonal([ones(1, 1), R₀(u)])

raw"Rotation matrix in $\mathbf{R}^3$, around Up (i.e., Z) axis."
R₃(u::Vector) = BlockDiagonal([R₀(u), ones(1, 1)])

"""Discrete sine wave"""
σ(n::Int) = sin(π * n / 2)

# --------------------- Control problem ---------------------

"""Shadow convex hull linear vector component"""
function a(panel::Panel, crop::Crop, sun::Sun, vertex_index)
    a₁ =
        0.5panel.width *
        panel.depth *
        [cos(sun.azimuth - panel.azimuth) / tan(sun.elevation), 1]
    a₂₁ =
        cos(sun.azimuth) / tan(sun.elevation) * (panel.pos[1] - crop.pos[1]) +
        sin(sun.azimuth) / tan(sun.elevation) * (crop.pos[2] - panel.pos[2])
    a₂₂ =
        cos(panel.azimuth) * (panel.pos[1] - crop.pos[1]) +
        sin(panel.azimuth) * (crop.pos[2] + panel.pos[2]) -
        panel.pos[3] * sin(sun.azimuth - panel.azimuth) / tan(sun.elevation)
    return a₁ + σ(vertex_index - 1) * panel.depth * [a₂₁, a₂₂]
end

"""Shadow convex hull linear scalar component"""
function b(panel::Panel, crop::Crop, sun::Sun, vertex_index)
    b₁ = panel.pos[3] * cos(sun.azimuth - panel.azimuth) / tan(sun.elevation)
    b₂ = sin(panel.azimuth) * (crop.pos[1] - panel.pos[1]) +
         cos(panel.azimuth) * (crop.pos[2] - panel.pos[2])
    return panel.width * σ(vertex_index) * (b₁ + b₂)
end

# If the returned expression is nonnegative, the crop is inside the 
# halfplane identified by the vertex_index
g(panel::Panel, crop::Crop, sun::Sun, vertex_index::Int, tilt_vec) =
    a(panel, crop, sun, vertex_index)' * tilt_vec + b(panel, crop, sun, vertex_index)

is_panel_shadowing_crop(panel::Panel, crop::Crop, sun::Sun, tilt_vec) =
    all(j -> g(panel, crop, sun, j, tilt_vec) >= 0, 1:4)

# Check (probabilistically) whether a panel can shadow a crop,
# given the sun position
function can_panel_shadow_crop(panel::Panel, crop::Crop, sun::Sun, num_samples::Int)
    # Deterministic pruning
    if (crop.pos .- panel.pos)' * light_beam(sun) <= 0
        return false
    else
    # Probabilistic pruning
        for angle in range(start = -π / 2, stop = π / 2, length = num_samples)
            if (
                is_panel_shadowing_crop(panel, crop, sun, u(angle)) &
                (proj_incidence_condition(panel, sun, u(angle)) >= 0)
            )
                return true
            end
        end
        return false
    end
end

# big-M and small-M for logical to integer constraint
# reformulation
function big_M(panel::Panel, crop::Crop, sun::Sun, vertex_index::Int)
    _a = a(panel, crop, sun, vertex_index)
    return -b(panel, crop, sun, vertex_index) + _a' * sign.(_a)
end

function small_M(panel::Panel, crop::Crop, sun::Sun, vertex_index::Int)
    _a = a(panel, crop, sun, vertex_index)
    return -b(panel, crop, sun, vertex_index) - _a' * sign.(_a)
end

# Optimal control open-loop
function control(sun::Sun, panels, crops, albedo, γ, num_samples=100, num_approx_points=6, ϵ=1e-4)
    print("[Sun @ $(sun.time)] Start constructing problem... ")
    prob = Model(SCIP.Optimizer)
    set_attribute(prob, "display/verblevel", 0)

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
    if γ > 0
      # Relaxed form
      @constraint(prob, tilt_vec[:, 1] .^ 2 .+ tilt_vec[:, 2] .^ 2 .<= 1)
    else
      # SOS2 formulation
      approx_range = range(start=-π / 2, stop=π / 2, length=num_approx_points) 
      approx_points = [sin.(approx_range) cos.(approx_range)]

      @variable(prob, 0 <= sos2_weight[1:length(panels), 1:num_approx_points] <= 1) # SOS2 weight
      @variable(prob, sos2_indicator[1:length(panels), 1:num_approx_points-1], Bin) # SOS2 indicator

      @constraints(
        prob,
        begin
          sum(sos2_weight, dims=2) .== 1
          sum(sos2_indicator, dims=2) .== 1
          sos2_weight[:,1:end-1] .+ sos2_weight[:,2:end] .>= sos2_indicator
          tilt_vec .== sos2_weight * approx_points
        end
      )
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
                    -g(panel, crop, sun, j, tilt_vec[i, :]) <=
                    big_M(panel, crop, sun, j) *
                    (1 - does_module_hp_shadow_crop[(i, j, k)])
                    -g(panel, crop, sun, j, tilt_vec[i, :]) >=
                    ϵ + (small_M(panel, crop, sun, j) - ϵ) * does_module_hp_shadow_crop[(i, j, k)]
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
        p(panel, sun, albedo, tilt_vec[i, :]) for (i, panel) in enumerate(panels)
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
    @objective(prob, Max, γ * total_power - (1 - γ) * reference_distance)
    optimize!(prob)
    println("Status: $(termination_status(prob)).")

    return value.(tilt_vec), value.(is_crop_shadowed)
end
