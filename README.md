# MOIA 
Under updating........
1. Introduction

 1.1 Dominating set (DSet)
 
Ore and Berge (1962) early introduced the concept of domination in graphs in which a subset of graph nodes covering all the nodes. Figure 1.1 shows an undirected graph in which  G= (V,E),  where V consists of vertices set and E consists of the edges set. A dominating set (shortly DS) is a subset of vertices D where  D⊆V such that for all u∈V-D, there exists a  v∈D, for which  uv∈E (we say that D dominates  V ). This means each vertex in V either belongs to the dominating set or it is adjacent to some vertices of the dominating set. Indeed, there are different variations of the dominating set.

1.2 	The MDSet problem

The MDSet problem is to find a DS of minimum size in a graph (Figures 1.4C and 1.4D). The MDSet represents a small-optimized subset in which any other node in the targeted graph must be at least adjacent to one node of the proposed MDSet (Nacher, J.C. and T. Akutsu,2012). Therefore, the MDSet covers all the nodes in the targeted network. The size of an MDSet in graph G is called the domination number of G and is denoted  γ(G). The MDSet problem is a fundamental problem in algorithmic graph theory. Moreover, MDSet problem is one of the central problems of combinatorial optimization that classified as NP-hard problem (Wuchty, S. 2014).

Nacher et al. (2014) presented the following ILP formulation to solve the MDSet problem.

![image](https://user-images.githubusercontent.com/53053110/115181576-8df5ac00-a0d8-11eb-908f-4bad2f919331.png)


1.3 The critical set

Nacher and Akutsu (2014) classified the nodes of graph depending on their belong to the MDSets generated by ILP-based model in to three types : critical (belongs to every MDSet), intermittent (may be belongs to one MDSet), and redundant (never belong to any MDSet). As showing in Figure 1.B, the node {3,6} represent the critical set that must be appear in every optimal MDSet (Figure 1.C and 1.D), {1,2,5,7,8,9,10} represent the intermittent nodes, and {4,11,12,14,15} represent the redundant nodes.



![image](https://user-images.githubusercontent.com/53053110/115179583-eeceb580-a0d3-11eb-9f72-ab7fef058709.png)


__________________________________________
2.Codes
__________________________________________
This Repository contains supplementary files for the project: “Constraint-Based Models for Dominating Protein Interaction Networks” 
The following Repository includes the following matlab functions:


2.1- MSKMDS_ILP : This function to solve graph minimun dominating set problem by ILP model using Mosek solver form matlab.

2.2-GRBMDS_ILP:This function to solve graph minimun dominating set problem by ILP model using Gurobi solver form matlab.

2.3-crds : This function to determine the critical , redundant  and intermittent  nodes using traditional method.

2.4-GetTwo_MDSets :This Function developed to solve BigMatrix by ILP solver to generate the most two different  minimum dominating sets connect matlab with mosek.

2.5-Get_MMDSets: This Function developed to solve BigMatrix by ILP solver to generate Multiple minimum dominating sets From the grpah adjacency Matrix.
    % The function also find the critical nodes as the intersection between the first generated Two MDSets.
 __________________________________________
3. DATASETS
 __________________________________________
    
3.1-Hint: Subfolder contains the protein Interaction networks (PPIN) dataset for Human and Yeast.
For Human PPI networks, we considered three different datasets obtained from H. Sapiens in the HINT database (version 3/10/2018) http://hint.yulab.org/download/. The first one of these datasets contains 63,684 high-quality binary protein (HHQBP) interactions between 12,815 human proteins. The second dataset contains 116,456 high-quality co-complex protein (HHQCP) interactions between 12,352 human proteins. However, a network of 180,140 combined protein (HCP) interactions between 15,744 human proteins is considered as the third dataset.

3.2-Bioplex : Subfolder contains the protein intwraction networks (PPIN) dataset for Human, Two versions of the protein interaction dataset of the BioPlex network [http://bioplex.hms.harvard.edu/)] were used. The first version, BIOPLEX1, had 23,744 proteins interactions between 7,637 proteins, and the second, BIOPLEX2, had 56,553 protein interactions between 10,883 proteins. Moreover, these two datasets with 80,297 protein interactions between 11,540 proteins were also combined as (BIOPLEX12). 

3.3-Liver : this subfolder includes 28,553 protein interactions between 7148 liver tissue proteins (LTP) collected in [X.-F. Zhang, L. Ou-Yang, D.-Q. Dai, M.-Y. Wu, Y. Zhu, and H. Yan, “Comparative analysis of housekeeping and tissue-specific driver nodes in human protein interaction networks,” BMC bioinformatics, vol. 17, no. 1, pp. 358, 2016.] was used.
