using Profile, BenchmarkTools, FlameGraphs

include("model.jl")

# Data paths
const PANELS_FILE = "data/simple-plant/modules.csv"
const CROPS_FILE = "data/simple-plant/crops.csv"
const SOLAR_FILE = "data/solar.csv"

# Optimization parameters
const GAMMA = 1
const ALBEDO = 0.01

panels, crops, suns, times = problem_data(PANELS_FILE, CROPS_FILE, SOLAR_FILE)

#= @btime control(suns[1], panels, crops, ALBEDO, GAMMA) =#
