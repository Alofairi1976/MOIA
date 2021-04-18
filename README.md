# MOIA 
Under updating........

__________________________________________
Codes
__________________________________________
This Repository contains supplementary files for the project: “Constraint-Based Models for Dominating Protein Interaction Networks” 
The following Repository includes the following matlab functions:


1- MSKMDS_ILP : This function to solve graph minimun dominating set problem by ILP model using Mosek solver form matlab.

2-GRBMDS_ILP:This function to solve graph minimun dominating set problem by ILP model using Gurobi solver form matlab.

3-crds : This function to determine the critical , redundant  and intermittent  nodes using traditional method.

4-GetTwo_MDSets :This Function developed to solve BigMatrix by ILP solver to generate the most two different  minimum dominating sets connect matlab with mosek.

5-Get_MMDSets: This Function developed to solve BigMatrix by ILP solver to generate Multiple minimum dominating sets From the grpah adjacency Matrix.
    % The function also find the critical nodes as the intersection between the first generated Two MDSets.
 __________________________________________
 DATASETS
 __________________________________________
    
6-Hint: Subfolder contains the protein intwraction networks (PPIN) dataset for Human and Yeast.
For Human PPI networks, we considered three different datasets obtained from H. Sapiens in the HINT database (version 3/10/2018) http://hint.yulab.org/download/. The first one of these datasets contains 63,684 high-quality binary protein (HHQBP) interactions between 12,815 human proteins. The second dataset contains 116,456 high-quality co-complex protein (HHQCP) interactions between 12,352 human proteins. However, a network of 180,140 combined protein (HCP) interactions between 15,744 human proteins is considered as the third dataset.

7-Bioplex : Subfolder contains the protein intwraction networks (PPIN) dataset for Human, Two versions of the protein interaction dataset of the BioPlex network [http://bioplex.hms.harvard.edu/)] were used. The first version, BIOPLEX1, had 23,744 proteins interactions between 7,637 proteins, and the second, BIOPLEX2, had 56,553 protein interactions between 10,883 proteins. Moreover, these two datasets with 80,297 protein interactions between 11,540 proteins were also combined as (BIOPLEX12). 

8-Liver : this subfolder includes 28,553 protein interactions between 7148 liver tissue proteins (LTP) collected in [X.-F. Zhang, L. Ou-Yang, D.-Q. Dai, M.-Y. Wu, Y. Zhu, and H. Yan, “Comparative analysis of housekeeping and tissue-specific driver nodes in human protein interaction networks,” BMC bioinformatics, vol. 17, no. 1, pp. 358, 2016.] was used.
