This directory contains the presence/absence tables from our analyses. The tables are in npz format which can be read in python. A simple script that would allow you to read the foldseek table using python and save it to a tab-delimited text file is provided below.

```python
import numpy as np
import scipy.sparse
import pandas as pd

#Matrices are coo (coordinate format) - that is they are sparse matrices
data = scipy.sparse.load_npz("Strep_foldseek_strain_by_gene.npz")
#One can see the size of the matrix by typing data.shape
print(data.shape)

#Let's read the labels
labels = []
with open ("Strep_foldseek_strain_by_gene.npz.labels.txt", 'r') as f:
    for line in f:
        labels.append(line.strip())

n_rows, n_cols = data.shape

#Convert the sparse matrix to array
m = data.toarray()
#Make it a data frame and add headers
df = pd.DataFrame(m, index=labels[:n_rows], columns=labels[n_rows:])

#See the first items
print(df.iloc[0:4,0:4])

#Print to a tab-delimited text file
df.to_csv('pangenomeMatrix.tsv', index=True, header=True, sep='\t')
```
