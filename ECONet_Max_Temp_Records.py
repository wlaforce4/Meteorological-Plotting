# Import libraries.
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as col
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import cartopy.io.shapereader as shprdr
import pygrib
from pyproj import CRS
import metpy

# Read in observed data.
obs = pd.read_csv(
    "TMAX_records.csv",
    names=[
        "station",
        "date",
        "currmax",
        "climomaxavg",
        "climomax",
        "prevmax",
        "monthmax",
    ],
)

# Calculate some new fields.
obs["tmaxanom"] = obs["currmax"] - obs["climomaxavg"]
obs["prevmaxdelta"] = obs["currmax"] - obs["prevmax"]
obs["monthmaxdelta"] = obs["currmax"] - obs["monthmax"]

# Read in metadata.
ECONet_meta = pd.read_csv("ECONet_metadata.csv")

# Merge our metadata to match that stations we have in the obs.
obs_meta = pd.merge(
    obs, ECONet_meta, how="left", left_on=["station"], right_on=["station"]
)

# Open the URMA grib file using the cfgrib engine
urma = xr.open_dataset("urma_nc.grb2", engine="cfgrib")

# We need to attach a coordinate reference system to our data for plotting.
# Quick way to get the correct coordinate reference system information from grib file.
grib = pygrib.open("urma_nc.grb2")
msg = grib.message(1)
crs_info = CRS(msg.projparams).to_cf()

# Assign both CRS info and a grid map to our data.This will allow us to index lat/lon points if desired.
urma = urma.metpy.assign_crs(
    {
        "semi_major_axis": crs_info.get("semi_major_axis"),
        "semi_minor_axis": crs_info.get("semi_minor_axis"),
        "grid_mapping_name": crs_info.get("grid_mapping_name"),
        "standard_parallel": [
            urma["t2m"].attrs["GRIB_Latin1InDegrees"],
            urma["t2m"].attrs["GRIB_Latin2InDegrees"],
        ],
        "latitude_of_projection_origin": urma["t2m"].attrs["GRIB_LaDInDegrees"],
        "longitude_of_central_meridian": urma["t2m"].attrs["GRIB_LoVInDegrees"],
    }
).metpy.assign_y_x()

# Creats CRS's for plotting.
# Create the CRS for the grib file
grid_crs = urma["t2m"].metpy.cartopy_crs

# Create the CRS for lat/lon points using the globe created from MetPy's 'assign_crs'
latlon_crs = ccrs.PlateCarree(globe=urma["t2m"].metpy.cartopy_globe)

# Create the CRS for the plotted map
map_crs = ccrs.AlbersEqualArea(
    central_longitude=-79.0,
    central_latitude=35.5,
    standard_parallels=(34.33, 36.16),
    globe=urma["t2m"].metpy.cartopy_globe,
)

# Since I'm looking for maxmimum temperatures in the URMA file, find the max value over the time dimension
urma_max = urma.max(dim="time")

# Create a map showing the gridded maximum temperatures from the URMA
# along with ECONet stations and their observed values.
# --------------------------------------------------------------------

fig = plt.figure(1, figsize=(11, 8.5))
ax = plt.subplot(1, 1, 1, projection=map_crs)
# Set the lat/lon extent to be plotted.
ax.set_extent([-84.275, -75.5, 33.2, 37.3])
# Set the background color.
ax.set_facecolor("whitesmoke")

# Read in shapefiles.
COUNTIES = cfeat.ShapelyFeature(shprdr.Reader("nc_cnty.shp").geometries(), latlon_crs)
STATE = cfeat.ShapelyFeature(shprdr.Reader("nc.shp").geometries(), latlon_crs)
NC_MASK = cfeat.ShapelyFeature(shprdr.Reader("nc_mask.shp").geometries(), latlon_crs)

# Plot the states and counties
ax.add_feature(NC_MASK, facecolor="whitesmoke", edgecolor="none", zorder=2)
ax.add_feature(COUNTIES, facecolor="none", edgecolor="grey", linewidth=0.15, zorder=3)
ax.add_feature(STATE, facecolor="none", edgecolor="black", linewidth=0.5, zorder=4)

# Plot the title and its location and font size
title_left = "Maximum Temperatures & ECONet Observations (°F) for Thursday Feb 09, 2023"
# title_right = "Valid: February 09, 2022"
plt.title(title_left, fontsize=11, fontweight="bold", loc="left")
# plt.title(title_right, fontsize=7, loc='right')

# Create a list of contour levels
clevs = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]

# Create the custom color map
cmap = col.ListedColormap(
    [
        "#4ac6fe",
        "#73d6fe",
        "#acfefe",
        "#30cec1",
        "#009895",
        "#125757",
        "#066d2c",
        "#31a254",
        "#74c376",
        "#a0d89a",
        "#d2febd",
        "#fefeb2",
        "#feec9f",
        "#fdd076",
        "#fdad2a",
        "#fc8c3c",
    ]
)

# Generates the colormap index based on discrete intervals
norm = col.BoundaryNorm(clevs, ncolors=cmap.N, extend="max")
# Set anything under our lowest clev value (i.e. 0.01 in this case) to be white
cmap.set_under(color="white")

# Plot the max temperatures as a contour fill
tmax = 1.8 * (urma_max["t2m"] - 273.15) + 32
s = ax.contourf(
    urma_max["longitude"],
    urma_max["latitude"],
    tmax,
    transform=latlon_crs,
    norm=norm,
    cmap=cmap,
    levels=clevs,
    zorder=1,
    extend="max",
)

# Get station value colors based on if they broke records or not.
val_colors = []
stat_count = 0
tot_stat = 0
for index, row in obs_meta.iterrows():
    if row["currmax"] > row["prevmax"]:
        clr = "red"
        stat_count += 1
        tot_stat += 1
    else:
        clr = "black"
        tot_stat += 1
    val_colors.append(clr)
val_cols = pd.DataFrame(val_colors, columns=["val_cols"])

# Label the ECONet stations with their respective value
values = []
for x, y, v, vx_o, vy_o, h_a, v_a, clr in zip(
    obs_meta["lon"],
    obs_meta["lat"],
    obs_meta["currmax"],
    obs_meta["vx_off"],
    obs_meta["vy_off"],
    obs_meta["h_a"],
    obs_meta["v_a"],
    val_cols["val_cols"],
):
    x_p, y_p = map_crs.transform_point(x, y, src_crs=latlon_crs)
    values.append(
        plt.text(
            x_p + vx_o,
            y_p + vy_o,
            v,
            ha=h_a,
            va=v_a,
            fontsize=6.5,
            path_effects=[pe.withStroke(linewidth=1.25, foreground="white")],
            fontweight="bold",
            color=clr,
            zorder=5,
        )
    )

# Plot the ECONet station points colored by the measured value
ax.scatter(
    obs_meta["lon"],
    obs_meta["lat"],
    c=obs_meta["currmax"],
    s=8,
    linewidths=0.5,
    edgecolor="black",
    transform=latlon_crs,
    zorder=5,
    norm=norm,
    cmap=cmap,
)

# Label the ECONet stations with their respective station names
stn_names = []
for x, y, stn, sy_o in zip(
    obs_meta["lon"], obs_meta["lat"], obs_meta["station"], obs_meta["sy_off"]
):
    x_p, y_p = map_crs.transform_point(x, y, src_crs=latlon_crs)
    stn_names.append(
        plt.text(
            x_p,
            y_p + sy_o,
            stn,
            va="center",
            ha="center",
            fontsize=3.75,
            color="black",
            zorder=5,
        )
    )

# Make new axis for the color bar and plot it on the bottom.
cax = fig.add_axes(
    [
        ax.get_position().x0 + 0.025,
        ax.get_position().y0 - 0.025,
        ax.get_position().width - 0.05,
        0.015,
    ]
)
cb = plt.colorbar(s, cax=cax, ticks=clevs, orientation="horizontal")
cb.ax.tick_params(labelsize=8)
cb.set_label("degrees F", size=8)

# Add logo to the plot.
logo = plt.imread("SCO_blue_335_v2.png")  # 150 dpi
fig.figimage(logo, 15, 101, zorder=6)  # 150 dpi

# Add some descriptive text to the plot.
ax.text(
    0.11,
    0.34,
    "Station values in red set a new record \n"
    + str(stat_count)
    + "/"
    + str(tot_stat)
    + "* stations broke their previous record \n*JRSP not included (<1 year of data)",
    fontsize=8,
    color="red",
    transform=ax.transAxes,
    va="bottom",
    ha="left",
)
ax.text(
    0.993,
    0.01,
    "Background contour fill may not represent actual conditions",
    fontsize=6,
    transform=ax.transAxes,
    va="bottom",
    ha="right",
)

# Save the image
plt.savefig(
    "Tmax_02092023.png",
    dpi=150,
    bbox_inches="tight",
    facecolor="white",
)
