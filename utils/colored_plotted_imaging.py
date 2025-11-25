import pandas as pd
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

DATA_PATH_TEMPLATE = 'D:/thesis_research/resources/{section_id}'
S10_SECTION_DIR = DATA_PATH_TEMPLATE.format(section_id='GSM8199188_ID61-ID62_S10')
S18_SECTION_DIR = DATA_PATH_TEMPLATE.format(section_id='GSM8199189_ID67-ID68_S18')
FIGURES_DIR = 'D:/thesis_research/resources/figures'

cells_s10 = pd.read_csv(f"{S10_SECTION_DIR}\\metadata_file.csv")




# Load data
section_dir = S10_SECTION_DIR
base_dir = Path(section_dir)
fov_csv = base_dir / "fov_positions_file.csv"
fov_df = pd.read_csv(fov_csv)

x_min = fov_df['x_global_px'].min()
y_min = fov_df['y_global_px'].min()
y_max = fov_df['y_global_px'].max()

# Load shapes
shapes_df = pd.read_csv("D:/thesis_research/resources/figures/Shapes_s10_bottom.csv")

shapes_df["napari_x"] = shapes_df["axis-1"]
shapes_df["napari_y"] = shapes_df["axis-0"]

# Transform cells
cells_s10["napari_x"] = cells_s10["CenterX_global_px"] - x_min
cells_s10["napari_y"] = y_max - cells_s10["CenterY_global_px"] + 4200

# Create polygons
region_polygons = {}
for shape_name, shape_data in shapes_df.groupby("index"):
    points = list(zip(shape_data["napari_x"], shape_data["napari_y"]))
    region_polygons[shape_name] = MplPath(points)


# Assign cells (vectorized)
def assign_all_regions(cells_df, polygons):
    points = list(zip(cells_df["napari_x"], cells_df["napari_y"]))
    regions = ["Unassigned"] * len(cells_df)

    for region_name, polygon in polygons.items():
        mask = polygon.contains_points(points)
        for i, is_inside in enumerate(mask):
            if is_inside and regions[i] == "Unassigned":
                regions[i] = region_name

    return regions


print("Assigning regions...")
cells_s10["region"] = assign_all_regions(cells_s10, region_polygons)

print("\nRegion counts:")
region_counts = cells_s10["region"].value_counts()
print(region_counts)

# Define colors
region_colors = {
    0: '#70AD47', 1: '#4472C4', 2: '#ED7D31',
    3: '#C55A5A', 'Unassigned': '#CCCCCC'
}

# Plot
fig, ax = plt.subplots(figsize=(20, 28), dpi=150)

# Plot each region
for region_name, color in region_colors.items():
    region_cells = cells_s10[cells_s10["region"] == region_name]
    if len(region_cells) > 0:
        ax.scatter(
            region_cells["napari_x"],
            region_cells["napari_y"],
            s=10,
            color=color,
            alpha=0.8,
            label=f"{region_name} ({len(region_cells):,})",
            edgecolors='none',
            rasterized=True
        )

# FIX: Set axis limits based on actual data
ax.set_xlim(cells_s10["napari_x"].min() - 1000, cells_s10["napari_x"].max() + 1000)
ax.set_ylim(cells_s10["napari_y"].min() - 1000, cells_s10["napari_y"].max() + 1000)

ax.set_aspect("equal")
ax.invert_yaxis()
ax.legend(loc='best', fontsize=14, markerscale=3)
ax.set_title("Cells Colored by Region", fontsize=18, weight='bold')
ax.set_xlabel("X (Napari coordinates)", fontsize=14)
ax.set_ylabel("Y (Napari coordinates)", fontsize=14)
ax.grid(True, alpha=0.2)

plt.tight_layout()
# plt.savefig("cells_colored_final.png", dpi=300, facecolor='white', bbox_inches='tight')
plt.show()

print("\n✓ Plot saved as 'cells_colored_final.png'")