# Spectral Permutation Recovery
R codes and data sets used in Ma, Cai and Li (2021, JASA)

Generally speaking, the spectral permutation recovery method takes a data set 
```math
SE = \frac{\sigma}{\sqrt{n}}
```
$X\in \R^{n\times p}$ with disordered columns as input. In the first step, each row is centred to have mean zero; in the second step, we apply SVD to the row-centred matrix and extract the first right singular vector $u_1$ associated to the largest singular value. The true order of columns are inferred based on the relative magnitude of the $u_1$ component.

The directory Supplement/Data/ contains the datasets used for Sections 6.2 and 6.3 of our main paper.
1. The files under the directory Subset12info contains the contig converges obtained from the 41 synthetic bacterial genomes. 
2. sample.csv contains the metadata of the samples used in Section 6.3.
3. all_PTR contains the estimated PTRs for all the samples at the 8th week obtained from DEMIC.


The directory Supplement/RCodes/ contains the R Codes used in this paper.
1. Wplot.R contains the R codes for Figure 2 in our online Supplementary Material.
2. Simulation_Synthetic.R contains the R codes for the simulations in Section 6.2 Figure 5 of our main paper.
3. Simulation_Supp.R contains the R codes for Figure 1 in our online Supplementary Material.
4. Simulation_01Loss.R contains the R codes for the simulations in Section 6.1 Table 1 of our main paper.
5. Simulation_Kendall.R contains the R codes for the simulations in Section 6.1 Figure 4 of our main paper.
