# LightGOEA
A light-weight python function for performing GO term enrichment analysis

# Documentation

```python
def GOEA(target_genes,GENE_SETS,goterms=None,fdr_thresh=0.25,p_thresh=1e-3): 
    """Performs GO term Enrichment Analysis using the hypergeometric distribution.
    
    Parameters
    ----------
    target_genes - array-like
        List of target genes from which to find enriched GO terms.
    GENE_SETS - dictionary
        Dictionary where the keys are GO terms and the values are lists of genes associated with each GO term.
        Ex: {'GO:0000001': ['GENE_A','GENE_B'],
             'GO:0000002': ['GENE_A','GENE_C','GENE_D']}
        Make sure to include all available genes that have GO terms in your dataset.
    goterms - array-list, optional, default None
        If provided, only these GO terms will be tested.
    fdr_thresh - float, optional, default 0.25
        Filter out GO terms with FDR q value greater than this threshold.
    p_thresh - float, optional, default 1e-3
        Filter out GO terms with p value greater than this threshold.
        
    Returns:
    -------
    enriched_goterms - pandas.DataFrame
        A Pandas DataFrame of enriched GO terms with FDR q values, p values, and associated genes provided.
    """
```    

# Usage
All you need is a target list of genes and a dictionary of GO terms with associated genes:
  Ex: 
  ```python
  GENE_SETS = {'GO:0000001': ['GENE_A','GENE_B'],
                   'GO:0000002': ['GENE_A','GENE_C','GENE_D']}
  ```
Make sure this dictionary contains all genes with available GO terms in your dataset.

To run GOEA, it is now simply:
```python
from light_goea import GOEA
enriched_goterms = GOEA(target_genes,GENE_SETS)
```
