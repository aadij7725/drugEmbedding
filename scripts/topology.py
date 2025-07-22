import json
import glob
import toponetx as tnx
from pathlib import Path
from typing import Set, List, Union

def load_json(path: Union[str, Path]):
    with open(path, 'r') as f:
        return json.load(f)

def build_drug_target_complex(
    drugs_path: Union[str, Path],
    targets_path: Union[str, Path],
    pairs_path: Union[str, Path],
    drug_group_glob: Union[str, Path] = "drug_groups_*.json",
    target_group_glob: Union[str, Path] = "target_groups_*.json",
) -> tnx.CombinatorialComplex:
    """
    Constructs a combinatorial complex from the provided drug/target jsons.

    Args:
        drugs_path: Path to drugs.json (list of drug IDs, as str or int)
        targets_path: Path to targets.json (list of target IDs, as str or int)
        pairs_path: Path to drug_target_pairs.json
        drug_group_glob: Glob pattern for drug_groups_*.json files
        target_group_glob: Glob pattern for target_groups_*.json files

    Returns:
        tnx.CombinatorialComplex object
    """
    # 0-cells: union of drugs and targets
    drug_ids = set(load_json(drugs_path))
    target_ids = set(load_json(targets_path))
    zero_cells = set(map(str, drug_ids)) | set(map(str, target_ids))  # ensure all are strings

    # 1-cells: all drug-target pairs
    pairs = load_json(pairs_path)
    one_cells = set()
    for pair in pairs:
        a = str(pair["drug_id"])
        b = str(pair["target_id"])
        one_cells.add(tuple(sorted([a, b])))

    # 2-cells: all groups from all group jsons
    two_cells = set()
    for group_path in glob.glob(str(drug_group_glob)):
        groups = load_json(group_path)
        for group in groups:
            # Each group is a 2-cell: the drug and its similar drugs
            cell = set([str(group["drug_id"])] + [str(d) for d in group["similar_drug_ids"]])
            if len(cell) > 1:
                two_cells.add(tuple(sorted(cell)))
    for group_path in glob.glob(str(target_group_glob)):
        groups = load_json(group_path)
        for group in groups:
            cell = set([str(group["target_id"])] + [str(t) for t in group["similar_target_ids"]])
            if len(cell) > 1:
                two_cells.add(tuple(sorted(cell)))

    # Build the combinatorial complex
    CC = tnx.CombinatorialComplex()
    for node in zero_cells:
        CC.add_cell(node, rank=0)
    for edge in one_cells:
        CC.add_cell(list(edge), rank=1)
    for face in two_cells:
        CC.add_cell(list(face), rank=2)
    return CC

# Example usage as a script
if __name__ == "__main__":
    drugs_path = "./outputs/drugs.json"
    targets_path = "./outputs/targets.json"
    pairs_path = "./outputs/drug_target_pairs.json"
    drug_group_glob = "./output/drug_groups_*.json"
    target_group_glob = "./outputs/target_groups_*.json"

    CC = build_drug_target_complex(
        drugs_path,
        targets_path,
        pairs_path,
        drug_group_glob,
        target_group_glob
    )
    print("Combinatorial complex built!")
    print(CC)