import numpy as np
import pandas as pd
import pvlib as pv

LAT, LON, ALT = 40.6902033, 16.7484484, 181 # Laterza 
START_DATE = "10-26-2024 00:00:00" # MM-DD-YYYY HH:MM:SS"
END_DATE = "10-27-2024 00:00:00"
FREQ = "15min"
TZ = "Europe/Rome"

location = pv.location.Location(LAT, LON, TZ, ALT)

time = pd.date_range(START_DATE, END_DATE, freq=FREQ)
data = pd.concat([location.get_clearsky(time), location.get_solarposition(time)], axis=1, join="inner")\
         .filter(["apparent_elevation", "azimuth", "dni", "dhi", "ghi"])\
         .query("(ghi > 0 | dni > 0 | dhi > 0) & apparent_elevation > 0")\
         .rename(columns={"apparent_elevation": "elevation"})

# Convert to radians and save
data[["azimuth", "elevation"]] = np.radians(data[["azimuth", "elevation"]])
data.to_csv("data/solar.csv", index_label="time")
