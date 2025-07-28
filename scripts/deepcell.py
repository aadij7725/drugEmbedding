import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import toponetx as tnx
import topology
import json
import glob
from pathlib import Path
import topoembedx as tex

# Create a cell complex object with a few cells
cc = topology.build_drug_target_complex(
    drugs_path="./sample/drugs.json",
    targets_path="./sample/targets.json",
    pairs_path="./sample/drug_target_pairs.json",
    drug_group_glob="./sample/drug_groups_*.json",
    target_group_glob="./sample/target_groups_*.json"
)

# Create a model
model = tex.DeepCell(dimensions=3)

# Fit the model to the combinatorial complex across all ranks

model.fit(cc, neighborhood_type="coadj", neighborhood_dim={"rank": 1, "via_rank": 0})

# Get the embeddings
embedded_points = model.get_embedding(get_dict=True)


# Prepare data for plotting
x = [embedded_points[cell][0] for cell in embedded_points]
y = [embedded_points[cell][1] for cell in embedded_points]
z = [embedded_points[cell][2] for cell in embedded_points]
cell_labels = [f"Cell {cell}" for cell in embedded_points]

# Plotting
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c="blue", label="Projected Points")

# Annotate the points to correspond with cells
#for i, label in enumerate(cell_labels):
#    ax.annotate(
#        label, (x[i], y[i], z[i]), textcoords="offset points", xyztext=(0, 10), ha="center"
#    )

# Label axes and add title
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.set_title("Projection of Cell Complex in 3D")

# Display the plot
plt.legend()
#plt.grid(True)
plt.show()