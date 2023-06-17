# Drug resitance mutation analysis in *Mycobacterium abscessus*

Data availability: 
- https://www.ncbi.nlm.nih.gov/bioproject/319839
- Davidson RM, Hasan NA, Epperson LE, Benoit JB, Kammlade SM, Levin AR, Calado de Moura V, Hunkins J, Weakly N, Beagle S, Sagel SD, Martiniano SL, Salfinger M, Daley CL, Nick JA, Strong M. 2021. Population Genomics of Mycobacterium abscessus from U.S. Cystic Fibrosis Care Centers. Ann Am Thorac Soc 18:1960-1969.

## 1) Notebooks 
- Jupyter notebook (Analysis_DR_mutations_MAB.ipynb) - contains the Python code used for the initial analyses. 
- R markdown (Analysis_DR_mutations.md) - contains the R code used for the statistical analyses and for visualizing phylogenetic trees.

## 2) Datasets
### Excel files:
- MAB_MIC_DR_mutations_DCC.xlsx - dataset with combined minimum inhibitory concentrations (MIC), drug resitance (DR), and dominant circulating clones (DCC) data for *M. abscessus* 
- MAB_dataset_clean.xlsx - clean dataset exported from initial analyses using Python 

### CSV files: 
- MAB_ATCC.csv - contains the label used for *M. abscessus* ATCC19977

Confusion matrices generated during the initial analyses  
- rrs_confusion_matrix.csv
- rrl_confusion_matrix.csv
- erm_confusion_matrix.csv
