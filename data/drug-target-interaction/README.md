### Steps to convert data from format of author (Dhanya Shridhar) to cli-interfaces format: 
make crd and ncrd seperate directories. All below steps apply to the 3 dirs of the 2 datasets

1. Copy all common files (validInteraction, ATCSimilarity, distSimilarity, seqSimilarity, ligandSimilarity, GOSimilarity, SideEffectSimilarity, chemicalSimilarity) into each split
2. All files is the directory is now eval partition.
3. For learn partition copy the next split into learn dir. Since this is module the 10th split gets 1st eval as its learn
4. rm all interacts_positives.csv and interacts_negatives.csv
5. Rename interacts and interactsids to truth and target respectively.
```
rename -n s/interacts.csv/interacts_truth.txt/g */*/*/interacts.csv
rename  s/interactsids.csv/interacts_target.txt/g */*/*/interactsids.csv
```

File names follow predicate_partitiontype.txt
