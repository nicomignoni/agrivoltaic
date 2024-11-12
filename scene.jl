using CSV, DataFrames, PyPlot

include("utils.jl")
using3D()

const COLOR_SHADOW_CROP = "black"
const COLOR_LIGHTING_CROP = "green"
const COLOR_PANEL = "blue"
const COLOR_SHADOW = "gray"
const COLOR_POLE = "black"

const SCALE_LIM = 2

const RESULTS_FILE = "2024-11-12T00:37:20.758"
const PANELS_FILE = "data/modules.csv"
const CROPS_FILE = "data/crops.csv"
const SOLAR_FILE = "data/solar.csv"

mkpath("media/$(RESULTS_FILE)")

# Get data
panels_data = CSV.read(PANELS_FILE, DataFrame)
results_data = CSV.read("results/$(RESULTS_FILE).csv", DataFrame)
crops_data = CSV.read(CROPS_FILE, DataFrame)

# Main plotting loop
for result in eachrow(results_data)
  
  # Plot base module name
  ax = subplot(projection="3d")
  ax.set_proj_type("persp")
  ax.view_init(elev=30, azim=230, roll=0)
  ax.set_axis_off()

  min_xlim, max_xlim = minimum(crops_data.pos_east), maximum(crops_data.pos_east)
  min_ylim, max_ylim = minimum(crops_data.pos_north), maximum(crops_data.pos_north)

  # Plot modules poles
  for i in [-1, 1]
    ax.stem(
      panels_data.pos_east + i * 0.5 * panels_data.width,
      panels_data.pos_north,
      panels_data.height,
      basefmt=" ",
      markerfmt=" ",
    )
  end

end 
