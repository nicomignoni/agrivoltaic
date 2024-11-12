import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import utils as u

# Settings
COLOR_SHADOW_CROP = "black"
COLOR_LIGHTING_CROP = "green"
COLOR_PANEL = "blue"
COLOR_SHADOW = "gray"
COLOR_POLE = "black"

SCALE_LIM = 2

RESULTS_FILE = "2024-11-12T17:45:24.203"
PANELS_FILE = "data/modules.csv"
CROPS_FILE = "data/crops.csv"
SOLAR_FILE = "data/solar.csv"

# Results directory
if not os.path.exists(f"media/{RESULTS_FILE}"):
    os.makedirs(f"media/{RESULTS_FILE}")

# Get data
panels_data = pd.read_csv(PANELS_FILE)
results_data = pd.read_csv(f"results/{RESULTS_FILE}.csv")
crops_data = pd.read_csv(CROPS_FILE)

# Base figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.set_proj_type("persp")
ax.view_init(elev=30, azim=230, roll=0)
ax.set_axis_off()

# Set axes limit
min_xlim, max_xlim = crops_data.pos_east.min(), crops_data.pos_east.max()
min_ylim, max_ylim = crops_data.pos_north.min(), crops_data.pos_north.max()

# Plot crops
crops_plot = ax.scatter(
    crops_data.pos_east,
    crops_data.pos_north,
    crops_data.height,
    depthshade=False
)

# Plot modules poles
for i in {-1, 1}:
   ax.stem(
       panels_data.pos_east + i * 0.5 * panels_data.width,
       panels_data.pos_north,
       panels_data.height,
       basefmt=" ",
       markerfmt=" "
   )

# Initialize panels and their shadow
panel_polygons = Poly3DCollection([np.zeros((4,3)) for i, _ in panels_data.iterrows()], color=COLOR_PANEL)
shadow_polygons = Poly3DCollection([np.zeros((4,3)) for i, _ in panels_data.iterrows()], color=COLOR_SHADOW)
ax.add_collection3d(panel_polygons)
ax.add_collection3d(shadow_polygons)

def update(t):
    result = results_data.loc[t,:]

    updated_panel_vertices = []
    updated_panel_shadow = []
    for i, panel in panels_data.iterrows():
        # Update vertices
        panel_vertices = u.vertices(
            panel.width,
            panel.depth,
            u.center(panel),
            u.u(panel.azimuth),
            u.u(result[f"planar_tilt_{i+1}"])
        )
        updated_panel_vertices.append(panel_vertices.T)

        # Update shadow 
        panel_shadow = u.projection_matrix(u.u(result.azimuth), u.u(result.elevation)) @ panel_vertices
        updated_panel_shadow.append(panel_shadow.T)

        panel_polygons.set_verts(updated_panel_vertices)
        shadow_polygons.set_verts(updated_panel_shadow)

        # Update crops color (depending on shading statdepending on shading statee)
        crops_color = [
            COLOR_SHADOW_CROP 
            if result[f"is_crop_{k+1}_shadowed"] == 1 
            else COLOR_LIGHTING_CROP 
            for k, _ in crops_data.iterrows()
        ]
        crops_plot.set_color(crops_color)

    return (panel_polygons, shadow_polygons, crops_plot)

anim = animation.FuncAnimation(fig, func=update, frames=results_data.shape[0])

ax.set_aspect("equal")
plt.show(block=False)
