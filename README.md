# MERIN
MERIN (mutations of epigenetic regulators in perturbed interactions) is an computational approach for identifying epigenetic regulator (ER) mutations associated with dysregulated
molecular interactions. With genomic profiles of tumor samples, molecular interaction data and interested epigenetic regulator genes, MERIN
can effectively identify the ERs potentially involved in the interactions.

## Usage
     MERIPT (exp.profile, Network, Mut, thr, pvalueCutoff, alpha,n.sim)
     
**Example**

``` r
load(file="./Project/IJob/MERIPT/exampleTest/exam_input.RData")

## Gene expression data
exp.profile[1:3,1:3]
#          TCGA-CG-4462 TCGA-CG-4306 TCGA-CG-4436
# AGRN         5.903019    5.2560569    4.3836984
# TNFRSF18     1.554199    1.5434400    0.6314070
# TNFRSF8      1.248957    0.2996354    0.1951491

## Molecular interaction data
Network[1:3,]
#    GeneA GeneB
#   ADAM10 EPHA3
#   ADAM23 ITGB3
#  ADCYAP1 GPR84

## Mutation data
Mut[1:3,1:3]
#              ACTR3B ACTR5 ATAD2
# TCGA-CG-4462      0     0     0
# TCGA-CG-4306      0     0     0
# TCGA-CG-4436      0     0     0

## Running the MERIPT
test_res <- MERIPT(exp.profile,Network,Mut)

## Showing the correlation results
test_res$MutPPI[1:3,]
#                     GeneA-B MutGene MutSamples DysSamples CooccuSamples    P
# MSTN;ACVR2B     MSTN;ACVR2B   ACTR5          3         40             2 0.04
# SCGB3A2;MARCO SCGB3A2;MARCO   ACTR5          3         24             2 0.00
# DUSP18;ITGA2   DUSP18;ITGA2   ATAD2         11          4             2 0.01
