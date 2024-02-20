# LouseEcomorphs
Title: **Parasite escape mechanisms drive morphological diversification in avian lice**
Authors: Stanislav Kolencik, Edward L. Stanley, Aswaj Punnath, Avery R. Grant, Jorge Dona, Kevin P. Johnson, and Julie M. Allen
2024, _Proceedings B_

Summary:
- nanoCT scan data were used to evaluate the morphological variability between 89 feather louse specimens from four ecomorph groups (body, generalist, head, and wing)
- specific morphological characters were quantified (head shape, proportional volume of mandibular muscularization, limb length)
- data were evaluated by statistical methods and geometric morphometric analysis, including phylogenomic data
- possible effects of louse-louse competition were explored
  
Outcomes: Feather lice have repeatedly evolved similar morphologies as a mechanism to escape host defenses and lice that co-occur with other louse genera on a host have shown greater morphological divergence.

Content:
1. **R script can be found in R markdown file "LouseEcomorph.md"**
2. **Files needed for these analyses are in .zip folder "LouseEcomorphsFiles.zip"**
   - .nts files contain 3D landmark information on head shape
   - Classifier.csv file is used in R analyses as a classifier file
   - various PCA scores (PaCA) .csv files are outputs from downstream analyses of landmark data
   - LouseTree_Binary.tre is a binary tree file -> some branches were split where our dataset includes landmarks for both female and male data (of the same specimen)
   - pruned_tree and tree_sorted are tree data from the downstream analyses of landmark data resulting from the inclusion of phylogenomic data
   - ClassifierCohab.csv and ClassifierLabel.csv files are useful for downstream analyses focusing on cohabiting vs. solo louse specimens or for including labels.
     
Zenodo: DOI 10.5281/zenodo.7557419
