using CSV, DataFrames, Plots, LaTeXStrings

include("model.jl")

# This is an addition
# ---------------- Read data ----------------

const GAMMA = 0.01
const ALBEDO = 0.1

const CROPS_FILE = "data/crops.csv"
const PANELS_FILE = "data/modules.csv"
const SOLAR_FILE = "data/solar.csv"

Plots.default(
  palette=:tab10,
  guidefontsize=8,
  xtickfontsize=6,            
  ytickfontsize=6,
  fontfamily="Computer Modern",
  legend=false,
  show=true
)

# rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
# rcParams["grid.alpha"] = 0.3
# rcParams["axes.spines.right"] = false
# rcParams["axes.spines.top"] = false
# rcParams["legend.frameon"] = false
# rcParams["ytick.labelsize"] = 6
# rcParams["xtick.labelsize"] = 6
# rcParams["font.size"] = 9
# rcParams["font.family"] = "sans"
# rcParams["font.sans-serif"] = ["Computer Modern Roman"]
# rcParams["text.usetex"] = true
# rcParams["text.latex.preamble"] = "\\usepackage{amsmath}";
# #     raw"\usepackage{amsfonts}", 
# #     raw"\usepackage{amssymb}",
# # ]

# ---------------- Read data ----------------

print("Loading data... ")
crops = [
  Crop([crop.pos_east, crop.pos_north, crop.height], crop.to_shadow) 
  for crop in CSV.read(CROPS_FILE, DataFrame) |> eachrow
]

panels = [
  Panel(
    panel.width, 
    panel.depth, 
    [panel.pos_east, panel.pos_north, panel.height],
    panel.azimuth
  )
  for panel in CSV.read(PANELS_FILE, DataFrame) |> eachrow
]

suns = [
  Sun(sun.time, sun.dni, sun.dhi, sun.ghi, sun.azimuth, sun.elevation)
  for sun in CSV.read(SOLAR_FILE, DataFrame) |> eachrow
]
println("Done")

# ---------------- Main loop ----------------

total_power = Vector{Real}(undef, length(suns))
light_coverage = Vector{Real}(undef, length(suns))
crops_to_light = sum(crop -> 1 - Int(crop.to_shadow), crops)
for (t, sun) in enumerate(suns)
  tilt_vec, is_crop_shadowed = control(sun, panels, crops, ALBEDO, GAMMA)

  total_power[t] = sum(
    panel_power(panel, sun, ALBEDO, tilt_vec[i,:])
    for (i, panel) in enumerate(panels)
  )
  light_coverage[t] = sum(1 .- is_crop_shadowed) / crops_to_light
end

# ---------------- Plots base ----------------

times = [sun.time for sun in suns]
# fig = subplots(figsize=(3.5, 2))
plot(1e-3total_power, size=(336, 200), xlabel="Time")
plot!(twinx(), light_coverage, color="tab:orange")
# ax.plot()
