using CSV, DataFrame, PyPlot

include("model.jl")

@pyimport mpl_toolkits
Poly3DCollection = mpl_toolkits.mplot3d.art3d.Poly3DCollection

# Colors
const COLOR_SHADOW_CROP = "black"
const COLOR_LIGHTING_CROP = "green"
const COLOR_PANEL = "blue"
const COLOR_SHADOW = "gray"
const COLOR_POLE = "black"

# Base figure
fig = figure()
ax = fig.add_subplot(projection="3d")
ax.set_proj_type("persp")
ax.view_init(elev=30, azim=230, roll=0)
ax.set_axis_off()

# Set axes limit
min_xlim, max_xlim = Inf, 0
min_ylim, max_ylim = Inf, 0
for crop in crops
    min_xlim, max_xlim = min(min_xlim, crop.pos[1]), max(max_xlim, crop.pos[1])
    min_ylim, max_ylim = min(min_ylim, crop.pos[2]), max(max_ylim, crop.pos[2])
end

# Plot crops
crops_pos = stack(crop -> crop.pos, crops)
crops_plot = ax.scatter(crops_pos[1,:], crops_pos[2,:], crops_pos[3,:], depthshade=false)

# Plot modules poles
panels_pos = stack(panel -> panel.pos, panels)
panels_width = stack(panel -> panel.width, panels)
for i in [-1, 1]
   ax.stem(
       panels_pos[1,:] + i * 0.5panels_width,
       panels_pos[2,:],
       panels_pos[3,:],
       basefmt=" ",
       markerfmt=" "
   )
end

# Initialize panels and their shadow
panel_polygons = Poly3DCollection([zeros(4,3) for _ in 1:length(panels)], color=COLOR_PANEL)
shadow_polygons = Poly3DCollection([zeros(4,3) for _ in 1:length(panels)], color=COLOR_SHADOW)
ax.add_collection3d(panel_polygons)
ax.add_collection3d(shadow_polygons)

function update_animation(t)
    sun = suns[t]

    updated_panel_vertices, updated_panel_shadow = [], []
    for panel in panels
        push!(
            updated_panel_vertices,
            
        )
    end
end
