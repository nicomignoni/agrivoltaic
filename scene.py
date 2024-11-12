import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import utils as u

# Settings
COLOR_SHADOW_CROP = "black"
COLOR_LIGHTING_CROP = "green"
COLOR_PANEL = "blue"
COLOR_SHADOW = "gray"
COLOR_POLE = "black"

SCALE_LIM = 2

RESULTS_FILE = "2024-11-12T10:40:22.187"
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

# Main plotting loop
for _, result in results_data.iterrows():

    # Plot base module
    ax = plt.figure().add_subplot(projection="3d")
    ax.set_proj_type("persp")
    ax.view_init(elev=30, azim=230, roll=0)
    ax.set_axis_off()

    min_xlim, max_xlim = crops_data.pos_east.min(), crops_data.pos_east.max()
    min_ylim, max_ylim = crops_data.pos_north.min(), crops_data.pos_north.max()

    # Plot modules poles
    for i in {-1, 1}:
        ax.stem(
            panels_data.pos_east + i * 0.5 * panels_data.width,
            panels_data.pos_north,
            panels_data.height,
            basefmt=" ",
            markerfmt=" ",
        )

    panels, shadows = [], []
    for i, panel in panels_data.iterrows():
        panel_center = np.array([[panel.pos_east, panel.pos_north, panel.height]]).T

        # Create panel's vertices
        panel_vertices = u.vertices(
            panel.width,
            panel.depth,
            panel_center,
            u.u(panel.azimuth),
            u.u(result[f"planar_tilt_{i+1}"]),
        )
        panels.append(panel_vertices.T)

        # Create panel's shadow
        panel_shadow = (
            u.projection_matrix(u.u(result.azimuth), u.u(result.elevation)) @ panel_vertices
        )
        shadows.append(panel_shadow.T)

    ax.add_collection3d(Poly3DCollection(panels, color=COLOR_PANEL))
    ax.add_collection3d(Poly3DCollection(shadows, color=COLOR_SHADOW))

    # Plot crops
    ax.scatter(
        crops_data.pos_east,
        crops_data.pos_north,
        crops_data.height,
        c=[COLOR_SHADOW_CROP if result[f"is_crop_{k+1}_shadowed"] == 1 else COLOR_LIGHTING_CROP for k, _ in crops_data.iterrows()],
        depthshade=False
    )

    ax.set_xlim(min_xlim - SCALE_LIM, max_xlim + SCALE_LIM)
    ax.set_ylim(min_ylim - SCALE_LIM, max_ylim + SCALE_LIM)
    ax.set_aspect("equal")
    plt.savefig(
        f"media/{RESULTS_FILE}/sun-at-{result.time.replace(" ", "T")}.png", dpi=300,
        bbox_inches="tight",
    )

plt.show(block=False)
