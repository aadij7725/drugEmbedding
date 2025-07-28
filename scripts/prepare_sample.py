import scipy.io
import numpy as np
import json
import os

def prompt_overwrite(path):
    if os.path.exists(path):
        response = input(f"File '{path}' already exists. Overwrite? [y/N]: ").strip().lower()
        return response == 'y'
    return True

# Load the .mat file
mat = scipy.io.loadmat('../data/drug-target-interaction/Perlman_Data.mat')

# Extract IDs
drug_ids = mat['ID_drugbank_drugs'].flatten()
target_ids = mat['ID_entrez_targets'].flatten()
interactions = mat['Interactions_Matrix']

# Compute 5% counts
n_drugs = len(drug_ids)
n_targets = len(target_ids)
n_drugs_5pct = max(1, n_drugs // 20)
n_targets_5pct = max(1, n_targets // 20)

drug_ids_5pct = drug_ids[:n_drugs_5pct]
target_ids_5pct = target_ids[:n_targets_5pct]

drug_ids_5pct_set = set(str(d) for d in drug_ids_5pct)
target_ids_5pct_set = set(str(t) for t in target_ids_5pct)

# 1. drugs.json
drugs_json_path = 'drugs.json'
if prompt_overwrite(drugs_json_path):
    with open(drugs_json_path, 'w') as f:
        json.dump([str(d) for d in drug_ids_5pct], f, indent=2)

# 2. targets.json
targets_json_path = 'targets.json'
if prompt_overwrite(targets_json_path):
    with open(targets_json_path, 'w') as f:
        json.dump([int(t) for t in target_ids_5pct], f, indent=2)

# 3. drug_target_pairs.json
pairs_json_path = 'drug_target_pairs.json'
pairs = []
for i, d_id in enumerate(drug_ids_5pct):
    for j, t_id in enumerate(target_ids_5pct):
        pairs.append({'drug_id': str(d_id), 'target_id': int(t_id), 'interaction': int(interactions[i, j])})
if prompt_overwrite(pairs_json_path):
    with open(pairs_json_path, 'w') as f:
        json.dump(pairs, f, indent=2)

# 4. Drug similarity matrices (groups only for first 5% drugs, similars restricted to first 5%)
drug_sims = {
    'ATCHier': mat['DrugSim_ATCHierDrugsCommonSimilarityMat'],
    'chemical': mat['DrugSim_chemicalDrugsCommonSimilarityMat'],
    'ligandJaccard': mat['DrugSim_ligandJaccardDrugsCommonSimilarityMat'],
    'newCMapJaccard': mat['DrugSim_newCMapJaccardDrugsCommonSimilarityMat'],
    'SideEffect': mat['DrugSim_pSideEffectDrugsCommonSimilarityMat']
}

drug_id_strs = [str(d) for d in drug_ids]
drug_id_to_idx = {str(d): i for i, d in enumerate(drug_ids)}

for sim_name, sim_mat in drug_sims.items():
    groups = []
    for i, d_id in enumerate(drug_ids_5pct):
        sim_indices = np.where((sim_mat[i] > 0) & (np.arange(n_drugs) != i))[0]
        similar_drugs = [str(drug_ids[idx]) for idx in sim_indices if str(drug_ids[idx]) in drug_ids_5pct_set]
        groups.append({'drug_id': str(d_id), 'similar_drug_ids': similar_drugs})
    out_path = f'drug_groups_{sim_name}.json'
    if prompt_overwrite(out_path):
        with open(out_path, 'w') as f:
            json.dump(groups, f, indent=2)

# 5. Target similarity matrices (groups only for first 5% targets, similars restricted to first 5%)
target_sims = {
    'dist': mat['TargetSim_distTargetsCommonSimilarityMat'],
    'GO': mat['TargetSim_GOTargetsCommonSimilarityMat'],
    'seq': mat['TargetSim_seqTargetsCommonSimilarityMat']
}

target_id_ints = [int(t) for t in target_ids]
target_id_to_idx = {int(t): i for i, t in enumerate(target_ids)}

for sim_name, sim_mat in target_sims.items():
    groups = []
    for i, t_id in enumerate(target_ids_5pct):
        sim_indices = np.where((sim_mat[i] > 0) & (np.arange(n_targets) != i))[0]
        similar_targets = [int(target_ids[idx]) for idx in sim_indices if str(target_ids[idx]) in target_ids_5pct_set]
        groups.append({'target_id': int(t_id), 'similar_target_ids': similar_targets})
    out_path = f'target_groups_{sim_name}.json'
    if prompt_overwrite(out_path):
        with open(out_path, 'w') as f:
            json.dump(groups, f, indent=2)

print("First 5% data files written (or skipped if not confirmed).")