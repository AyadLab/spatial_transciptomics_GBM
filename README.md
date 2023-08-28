# spatial_transciptomics_GBM

DATA INFO :

- We used one of GEO dataset,where they collected five DMG patients, five GBM patients (includes two secondary GBM), and one peritumor samples (total 11 patient samples) and they performed short-read spatial transcriptomic sequencing on the tissue sections using the 10x Visium platform.
- link to the data is given below
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194329)


DESCRIPTION :

- Our main focus is to investigate Casein Kinase 2 (CK2), which plays an active role in many cancers by controlling various signaling pathways. A lot of clinical trials and studies on mouse models observed that CK2 inhibition is needed in tumor cells as they use CK2 for survival. We want to investigate whether CK2 is expressed or co-expressed with certain other genes and also see if there are any implications of this protein in GBM tumors.Once we determine the spatial expression of CK2, we could see co-expressed genes,by checking whether other genes that are implicated in GBM are expressed in similar patterns.
- After this,we applied deconvolution on our integrated seurat object which consists of information on multiple assays like spatial,SCT,RNA and integrated.With deconvolutuon,we try to see what cell types at what proportion are present in our single RNA-seq data
- With Enrichment analysis,we try to achieve the relation between genes,clusters and their pathways to see any biological significance with it.


# packages used :
1. Seurat : to deal with spatial images.
2. STdeconvolve (unsupervised algorithm) : deconvolution of images(integrated)
3. msigdbr : database(human) used for enrichment analysis to find spatail pathways of each brain tissue using GESECA analysis
