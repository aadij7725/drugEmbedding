import scipy.io
import numpy as np
import pandas as pd
import json

# Load the .mat file
mat = scipy.io.loadmat('../data/drug-target-interaction/Perlman_Data.mat')

# Extract IDs
drug_ids = mat['ID_drugbank_drugs'].flatten()
target_ids = mat['ID_entrez_targets'].flatten()
interactions = mat['Interactions_Matrix']

# 1. Save drugs.json and targets.json
with open('drugs.json', 'w') as f:
    json.dump([str(d) for d in drug_ids], f, indent=2)

with open('targets.json', 'w') as f:
    json.dump([int(t) for t in target_ids], f, indent=2)

# 2. Save all drug-target pairs with interaction label
pairs = []
for i, d_id in enumerate(drug_ids):
    for j, t_id in enumerate(target_ids):
        pairs.append({'drug_id': str(d_id), 'target_id': int(t_id), 'interaction': int(interactions[i, j])})
with open('drug_target_pairs.json', 'w') as f:
    json.dump(pairs, f, indent=2)

# 3. Grouping by similarity matrices

# Drug similarity matrices
drug_sims = {
    'ATCHier': mat['DrugSim_ATCHierDrugsCommonSimilarityMat'],
    'chemical': mat['DrugSim_chemicalDrugsCommonSimilarityMat'],
    'ligandJaccard': mat['DrugSim_ligandJaccardDrugsCommonSimilarityMat'],
    'newCMapJaccard': mat['DrugSim_newCMapJaccardDrugsCommonSimilarityMat'],
    'SideEffect': mat['DrugSim_pSideEffectDrugsCommonSimilarityMat']
}

for sim_name, sim_mat in drug_sims.items():
    groups = []
    for i, d_id in enumerate(drug_ids):
        # Exclude self-similarity, get all drugs with similarity > 0
        sim_indices = np.where((sim_mat[i] > 0) & (np.arange(len(drug_ids)) != i))[0]
        similar_drugs = [str(drug_ids[idx]) for idx in sim_indices]
        groups.append({'drug_id': str(d_id), 'similar_drug_ids': similar_drugs})
    with open(f'drug_groups_{sim_name}.json', 'w') as f:
        json.dump(groups, f, indent=2)

# Target similarity matrices
target_sims = {
    'dist': mat['TargetSim_distTargetsCommonSimilarityMat'],
    'GO': mat['TargetSim_GOTargetsCommonSimilarityMat'],
    'seq': mat['TargetSim_seqTargetsCommonSimilarityMat']
}

for sim_name, sim_mat in target_sims.items():
    groups = []
    for i, t_id in enumerate(target_ids):
        # Exclude self-similarity, get all targets with similarity > 0
        sim_indices = np.where((sim_mat[i] > 0) & (np.arange(len(target_ids)) != i))[0]
        similar_targets = [int(target_ids[idx]) for idx in sim_indices]
        groups.append({'target_id': int(t_id), 'similar_target_ids': similar_targets})
    with open(f'target_groups_{sim_name}.json', 'w') as f:
        json.dump(groups, f, indent=2)

print("All JSON files written.")