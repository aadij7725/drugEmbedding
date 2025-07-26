import matplotlib.pyplot as plt
import toponetx as tnx
import topology
import json
import glob
from pathlib import Path
import topoembedx as tex

# Create a cell complex object with a few cells
cc = topology.build_drug_target_complex(
    drugs_path="./outputs/drugs.json",
    targets_path="./outputs/targets.json",
    pairs_path="./outputs/drug_target_pairs.json",
    drug_group_glob="./outputs/drug_groups_*.json",
    target_group_glob="./outputs/target_groups_*.json"
)

# Create a model
model = tex.Cell2Vec(dimensions=2)

# Fit the model to the combinatorial complex across all ranks

model.fit(cc, neighborhood_type="adj", neighborhood_dim={"rank": -1, "via_rank": -1})

# note that "via_rank" is ignored here and only considered when the complex is Combintatorial complex

# Get the embeddings
embedded_points = model.get_embedding(get_dict=True)


# Prepare data for plotting
x = [embedded_points[cell][0] for cell in embedded_points]
y = [embedded_points[cell][1] for cell in embedded_points]
cell_labels = [f"Cell {cell}" for cell in embedded_points]

# Plotting
plt.figure(figsize=(10, 8))
plt.scatter(x, y, c="blue", label="Projected Points")

# Annotate the points to correspond with cells
for i, label in enumerate(cell_labels):
    plt.annotate(
        label, (x[i], y[i]), textcoords="offset points", xytext=(0, 10), ha="center"
    )

# Label axes and add title
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Projection of Cell Complex in 2D")

# Display the plot
plt.legend()
plt.grid(True)
plt.show()