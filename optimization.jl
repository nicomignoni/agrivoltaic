using Dates, LinearAlgebra, CSV, DataFrames, JuMP, MosekTools

include("utils.jl")

# Settings
const S = 10
const EPSILON = 1e-4
const GAMMA = 0.01
const ALBEDO = 0.1

const CROPS_FILE = "data/crops.csv"
const PANELS_FILE = "data/modules.csv"
const SOLAR_FILE = "data/solar.csv"

# Read data 
print("Loading data... ")
crops_data = CSV.read(CROPS_FILE, DataFrame);
panels_data = CSV.read(PANELS_FILE, DataFrame);
solar_data = CSV.read(SOLAR_FILE, DataFrame);
println("Done")

# The results dataframe collects sun related data, the 
# tilt for each module, and the shadowing status of
# each crop
results = hcat(
    solar_data,
    DataFrame(["planar_tilt_$(i)" => zeros(nrow(solar_data)) for i = 1:nrow(panels_data)]),
    DataFrame(["is_crop_$(k)_shadowed" => zeros(nrow(solar_data)) for k = 1:nrow(crops_data)])
 )

# Main
for (t, sun) in eachrow(solar_data) |> enumerate
    print("[Sun @ $(sun.time)] Start constructing problem... ")
    prob = Model(MosekTools.Optimizer)
    set_attribute(prob, "QUIET", true)
    
    # Consider a pair (panel, crop) only if the panel can actually
    # shadow the crop
    I = [
      (i, k) for (i, panel) in eachrow(panels_data) |> enumerate for
      (k, crop) in eachrow(crops_data) |> enumerate if
      can_panel_shadow_crop(panel, crop, sun, S)
     ]

    # Planar tilt vector, with bounds
    @variable(prob, z - 2 <= tilt_vec[1:nrow(panels_data), z = 1:2] <= 1)

    # Is the crop position shadowed?
    @variable(prob, is_crop_shadowed[1:nrow(crops_data)], Bin, start=0)

    # Fundamental trigonometric (relaxed) equality constraint
    if GAMMA > 0
        @constraint(prob, tilt_vec[:, 1] .^ 2 .+ tilt_vec[:, 2] .^ 2 .<= 1)
    end

    # Avoid light to hit the back of the panel 
    # (negative irradiance flux -> negative power) 
    @constraint(
        prob,
        [(i, panel) in eachrow(panels_data) |> enumerate],
        proj_incidence_condition(panel, sun, tilt_vec[i,:]) >= 0
    )

    num_modules_shadowing_crop = zeros(AffExpr, nrow(crops_data))

    # Does the panel shadow the crop position?
    does_module_shadow_crop = Dict{NTuple{2,Int}, VariableRef}()

    # Does the panel's halfplane shadows the crop position?
    does_module_hp_shadow_crop = Dict{NTuple{3,Int}, VariableRef}()
    for (i, k) in I
        does_module_shadow_crop[(i, k)] = @variable(prob, binary=true)

        panel = panels_data[i, :]
        crop = crops_data[k, :]

        # Shadow indicators definitions' constraints
        for j = 1:4
            does_module_hp_shadow_crop[(i, j, k)] = @variable(prob, binary=true)
            @constraints(
                prob,
                begin
                  - shadowing_condition(panel, crop, sun, j, tilt_vec[i,:]) <=
                    big_M(panel, crop, sun, j) * (1 - does_module_hp_shadow_crop[(i, j, k)])
                  - shadowing_condition(panel, crop, sun, j, tilt_vec[i,:]) >=
                    EPSILON + (small_M(panel, crop, sun, j) - EPSILON) * does_module_hp_shadow_crop[(i, j, k)]
                end
            )
        end

        num_modules_shadowing_crop[k] += does_module_shadow_crop[(i, k)]

        # Is the crop contained in all the modules' halfplanes (i.e., is 
        # the module shadowing the crop)? 
        @constraint(
            prob,
            sum(j -> does_module_hp_shadow_crop[(i, j, k)], 1:4) <=
            does_module_shadow_crop[(i, k)] + 3
        )
        @constraint(
            prob,
            sum(j -> does_module_hp_shadow_crop[(i, j, k)], 1:4) >=
            4does_module_shadow_crop[(i, k)]
        )
    end

    # Does any module shadow the crop?
    @constraints(
        prob,
        begin
            num_modules_shadowing_crop .<= nrow(panels_data) * is_crop_shadowed
            num_modules_shadowing_crop .>= is_crop_shadowed
        end
    )

    # Objective function terms
    total_power = sum(
        panel_power(panel, sun, ALBEDO, tilt_vec[i, :]) for
        (i, panel) in eachrow(panels_data) |> enumerate
    )

    # In JuMP, the l1 norm is explicitly defined throught the canonical cone formulation
    shadow_reference = stack(crop.to_shadow for crop in eachrow(crops_data))
    @variable(prob, reference_distance)
    @constraint(
        prob,
        [reference_distance; is_crop_shadowed - shadow_reference] in
        MOI.NormOneCone(1 + nrow(crops_data))
    )
    print("Done. ")

    # Solve
    @objective(prob, Max, GAMMA * total_power - (1 - GAMMA) * reference_distance)
    optimize!(prob)
    println("Status: $(termination_status(prob)).")

    # Update results
    for (i, vec) in eachrow(tilt_vec) |> enumerate
      results[t, "planar_tilt_$(i)"] = angle(vec) |> value
    end
    for (k, status) in is_crop_shadowed |> enumerate
      results[t, "is_crop_$(k)_shadowed"] = round(status |> value)
    end
end

# Save results
CSV.write("results/$(Dates.now()).csv", results)
