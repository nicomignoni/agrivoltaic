using CSV, DataFrames, LaTeXStrings, PyCall, PyPlot, ThreadsX

include("model.jl")

# Optimization parameters
const GAMMA = 1
const ALBEDO = 0.01

# Data paths
const CROPS_FILE = "data/la-svolta-plant/crops.csv"
const PANELS_FILE = "data/la-svolta-plant/modules.csv"
const SOLAR_FILE = "data/solar.csv"

# PyPlot settings
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["figure.constrained_layout.use"] = true
rcParams["grid.alpha"] = 0.3
# rcParams["axes.spines.right"] = false
rcParams["axes.spines.top"] = false
rcParams["legend.frameon"] = false
rcParams["ytick.labelsize"] = 6
rcParams["xtick.labelsize"] = 6
rcParams["font.size"] = 8
rcParams["font.family"] = "sans"
rcParams["font.sans-serif"] = ["Computer Modern Roman"]
rcParams["text.usetex"] = true
rcParams["text.latex.preamble"] = "\\usepackage{amsmath}";
#     raw"\usepackage{amsfonts}", 
#     raw"\usepackage{amssymb}",
# ]

formatter = matplotlib.dates.DateFormatter("%H:%M")
locator = matplotlib.dates.HourLocator(interval=2)

# Read data
panels, crops, suns, times = problem_data(PANELS_FILE, CROPS_FILE, SOLAR_FILE)
optimal_tilt_vecs = Array{Real}(undef, (length(suns), length(panels), 2))

optimal_total_power = Vector{Real}(undef, length(suns))
optimal_light_coverage = Vector{Real}(undef, length(suns))

vertical_policy_total_power = Vector{Real}(undef, length(suns))
vertical_policy_light_coverage = Vector{Real}(undef, length(suns))

horizontal_policy_total_power = Vector{Real}(undef, length(suns))
horizontal_policy_light_coverage = Vector{Real}(undef, length(suns))

classical_policy_total_power = Vector{Real}(undef, length(suns))
classical_policy_light_coverage = Vector{Real}(undef, length(suns))

crops_reference = [crop.to_shadow for crop in crops]
#= ThreadsX.foreach(enumerate(suns)) do (t, sun) =#
for (t, sun) in enumerate(suns)
    # Optimal control policy
    optimal_tilt_vecs[t,:,:], is_crop_shadowed = control(sun, panels, crops, ALBEDO, GAMMA)
    optimal_total_power[t] = sum(panel_power(panel, sun, ALBEDO, optimal_tilt_vecs[t,i,:]) for (i, panel) in enumerate(panels))
    optimal_light_coverage[t] = (is_crop_shadowed .- crops_reference) .|> abs |> sum

    # 180-deg tilt
    vertical_tilt_vecs = [u(π/2) for  _ in panels]
    vertical_policy_total_power[t] = sum(panel_power(panel, sun, ALBEDO, vertical_tilt_vecs[i]) for (i,panel) in enumerate(panels))
    vertical_policy_light_coverage[t] = (
        [f_is_crop_shadowed(panels, crop, sun, vertical_tilt_vecs) for crop in crops] .- crops_reference
    ) .|> abs |> sum

    # 30-deg tilt
    classical_tilt_vecs = [u(π/6) for _ in panels]
    classical_policy_total_power[t] = sum(panel_power(panel, sun, ALBEDO, classical_tilt_vecs[i]) for (i,panel) in enumerate(panels))
    classical_policy_light_coverage[t] = (
        [f_is_crop_shadowed(panels, crop, sun, classical_tilt_vecs) for crop in crops] .- crops_reference
    ) .|> abs |> sum

    # 0-deg tilt
    horizontal_tilt_vecs = [u(0) for _ in panels]
    horizontal_policy_total_power[t] = sum(panel_power(panel, sun, ALBEDO, horizontal_tilt_vecs[i]) for (i,panel) in enumerate(panels))
    horizontal_policy_light_coverage[t] = (
        [f_is_crop_shadowed(panels, crop, sun, horizontal_tilt_vecs) for crop in crops] .- crops_reference
    ) .|> abs |> sum
end

# Optimal control plotting
fig, axs = subplots(2, figsize=(3.5, 2.5), sharex=true, height_ratios=[1, 1])

ax_right, ax_tilt = axs
ax_right.xaxis.set_major_formatter(formatter)
ax_right.xaxis.set_major_locator(locator)

ax_left = ax_right.twinx()
ax_right.plot(times, 1e-3optimal_total_power, color="tab:orange", label=L"\sum_{i \in \mathcal{M}} p_i")
ax_left.plot(
    times, optimal_light_coverage / length(crops), color="tab:green", 
    label=L"\frac{1}{|\mathcal{P}|}|\boldsymbol{\xi} - \boldsymbol{\xi}^\bullet|"
)

# ax_tilt.xaxis.set_major_formatter(ticker)
ax_tilt.set_xlabel("Time")
ax_tilt.set_ylabel("[Rad.]")
ax_tilt.set_xlim(times[1], times[end])
ax_tilt.set_ylim(-π/2 - 1, π/2 + 1)
ax_tilt.grid(true)

ax_right.set_ylabel("[kW]")
ax_right.grid(true)

ax_left.set_ylim(0, 1)

ax_left.legend(
    bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower right",
    borderaxespad=0, ncol=1
)
ax_right.legend(
    bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
    borderaxespad=0, ncol=1
)

tilt_points = stack(
    (time, angle(optimal_tilt_vecs[t,i,:]))
    for (t, time) in enumerate(times), (i, _) in enumerate(panels);
    dims=1
)
ax_tilt.scatter(tilt_points[:,1], tilt_points[:,2], s=0.5)
ax_tilt.fill_between(times, -π/2, π/2, color="gray", alpha=0.2)
ax_tilt.axhline(π/2, color="gray", ls="--")
ax_tilt.axhline(-π/2, color="gray", ls="--")

savefig("media/gamma-$(GAMMA).pdf")

# Comparison policies
fig, ax_right = subplots(figsize=(3.5, 1.8))

ax_right.xaxis.set_major_formatter(formatter)
ax_right.xaxis.set_major_locator(locator)

ax_left = ax_right.twinx()

power_curves = [
    (vertical_policy_total_power, "--", L"\text{Vertical modules policy, $\beta_i = \frac{\pi}{2}, \ \forall i \in \mathcal{M}$}"),
    (horizontal_policy_total_power, ":", L"\text{Horizontal modules policy, $\beta_i = 0, \ \forall i \in \mathcal{M}$"),
    (classical_policy_total_power, "-.", L"\text{Rule-of-thumb policy, $\beta_i = \frac{\pi}{6}, \ \forall i \in \mathcal{M}$}")
]

light_curves = [
    (vertical_policy_light_coverage, "--"),
    (horizontal_policy_light_coverage, ":"),
    (classical_policy_light_coverage, "-.")
]

for (curve, ls, label) in power_curves
    ax_right.plot(times, 1e-3curve, color="gray", ls=ls, label=label)
    ax_right.plot(times, 1e-3curve, color="tab:orange", ls=ls)
end

for (curve, ls) in light_curves
    ax_left.plot(times, curve / length(crops), color="tab:green", ls=ls)
end

ax_left.set_ylim(0, 1)

ax_right.set_xlabel("Time")
ax_right.set_ylabel("[kW]")

ax_right.legend(
    bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower right",
    borderaxespad=0, mode="expand", ncol=1, fontsize=6
)

ax_right.grid(true)
savefig("media/fix-tilt.pdf")

# Considered panel-crop pairs
spared_variable_ratio = Vector{Real}(undef, length(suns))
for (t, sun) in enumerate(suns)
    considered_panel_crop_pairs = sum(
        Int(can_panel_shadow_crop(panel, crop, sun, 100))
        for panel in panels, crop in crops
    )
    spared_variable_ratio[t] = considered_panel_crop_pairs / (length(crops) * length(panels))  
end

fig, ax = subplots(figsize=(3.5, 1.2))
ax.bar(times, 1e2spared_variable_ratio, linewidth=0, width=8e-3)
max_spared_variable_ratio = 1e2maximum(spared_variable_ratio)
ax.axhline(max_spared_variable_ratio, ls="--", lw=0.8, color="gray")
ax.text(times[2], 0.5 + max_spared_variable_ratio, "max. $(round(max_spared_variable_ratio))\\%", fontsize=6)

ax.xaxis.set_major_formatter(formatter)
ax.xaxis.set_major_locator(locator)

ax.set_xlabel("Time")
ax.set_ylabel(L"\frac{|\mathcal{I}|}{|\mathcal{M}||\mathcal{P}|} \ [\%]")

ax.set_ylim(0, 20)
ax.grid(true)

savefig("media/panel-crop-pairs.pdf")
