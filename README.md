# ITH Programs
### Program Description
This program implements the method described by Gavish et al.[6] to find the Intra-tumoral Heterogeneity (ITH) expression programs. The initial objects were a log normalized cell expression matrices with rows as genes and cells as columns.

1. For each one of our patient's cell expression matrices, we ran 7 different NMF
by using a k ranging from 5 to 11.
2. Based on the NMF score, each NMF Program was then summarized by its top
50 genes.
3. After multiple filters, the remaining NMF Programs were refered to as robust NMF Programs and then clustered in Meta-Programs (MPs).

### Clustering method
The clustering method is described in the figure: 
<p align="center">
<img src="clustering.png" alt="1" width="400"/>
</p> 




