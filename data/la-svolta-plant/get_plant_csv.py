import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plant_df = pd.read_csv("plant.csv")

# Crop data
crops_df = plant_df.query("width == depth")\
                   .drop(["width", "depth"], axis=1)
crops_df["to_shadow"] = True
crops_df["height"] = 0
crops_df["pos_north"] = -crops_df["pos_north"]

# Panel data
panels_df = plant_df.query("~(width == depth)")
panels_df["height"] = 2.5
panels_df["azimuth"] = np.pi
panels_df["pos_north"] = -panels_df["pos_north"]
panels_df["pos_east"] += 40

# Save as .csv
crops_df.to_csv("crops.csv")
panels_df.to_csv("modules.csv")

# Checking the overall plant
fig, ax = plt.subplots()
ax.scatter(crops_df["pos_east"], crops_df["pos_north"], s=1, color="tab:green")

panels_poly = []
for _, panel in panels_df.iterrows():
    verts = np.array([
        [panel.pos_east - 0.5*panel.width, panel.pos_north - 0.5*panel.depth],
        [panel.pos_east + 0.5*panel.width, panel.pos_north - 0.5*panel.depth],
        [panel.pos_east + 0.5*panel.width, panel.pos_north + 0.5*panel.depth],
        [panel.pos_east - 0.5*panel.width, panel.pos_north + 0.5*panel.depth],
    ])
    panels_poly.append(verts)
ax.add_collection(PolyCollection(panels_poly))

plt.autoscale(enable=True, axis='both')
ax.set_aspect("equal")
plt.show(block=False)

