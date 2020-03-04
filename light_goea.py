import numpy as np
import pandas as pd
from scipy.stats import rankdata

def log_factorial(n):
    return np.log(np.arange(1,n+1)).sum()
def log_binomial(n,k):
    return log_factorial(n) - (log_factorial(k) + log_factorial(n-k))

def GOEA(target_genes,GENE_SETS,df_key='GO',goterms=None,fdr_thresh=0.25,p_thresh=1e-3): 
    """Performs GO term Enrichment Analysis using the hypergeometric distribution.
    
    Parameters
    ----------
    target_genes - array-like
        List of target genes from which to find enriched GO terms.
    GENE_SETS - dictionary or pandas.DataFrame
        Dictionary where the keys are GO terms and the values are lists of genes associated with each GO term.
        Ex: {'GO:0000001': ['GENE_A','GENE_B'],
             'GO:0000002': ['GENE_A','GENE_C','GENE_D']}
        Make sure to include all available genes that have GO terms in your dataset.
        
        ---OR---
        
        Pandas DataFrame with genes as the index and GO terms values.
        Ex: 'GENE_A','GO:0000001',
            'GENE_A','GO:0000002',
            'GENE_B','GO:0000001',
            'GENE_B','GO:0000004',
            ...
        If `GENE_SETS` is a pandas DataFrame, the `df_key` parameter should be the name of the column in which
        the GO terms are stored.       
    df_key - str, optional, default 'GO'
        The name of the column in which GO terms are stored. Only used if `GENE_SETS` is a DataFrame.
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
    
    # identify all genes found in `GENE_SETS`
    all_genes = np.unique(np.concatenate(list(GENE_SETS.values())))
    all_genes = np.array(all_genes)
    
    if isinstance(GENE_SETS,pd.DataFrame):
        print('Converting DataFrame into dictionary')
        genes = GENE_SETS.index
        agt = GENE_SETS[df_key].values
        goterms = np.unique(agt)
        gs={}
        for go in goterms:
            gs[go] = np.unique(np.array(list(GENE_SETS.index[agt==go])))
        GENE_SETS = gs
    
    # if goterms is None, use all the goterms found in `GENE_SETS`
    if goterms is None:
        goterms = np.unique(list(GENE_SETS.keys()))
    else:
        goterms = goterms[np.in1d(goterms,np.unique(list(GENE_SETS.keys())))]
    
    # ensure that target genes are all present in `all_genes`
    _,ix = np.unique(target_genes,return_index=True)
    target_genes=target_genes[np.sort(ix)]
    target_genes = target_genes[np.in1d(target_genes,all_genes)]
    
    # N -- total number of genes
    N = all_genes.size

    probs=[]
    probs_genes=[]
    counter=0
    # for each go term,
    for goterm in goterms:
        if counter%1000==0:
            print(counter)
        counter+=1
        
        # identify genes associated with this go term
        gene_set = np.array(GENE_SETS[goterm])
        
        # B -- number of genes associated with this go term
        B = gene_set.size
        
        # b -- number of genes in target associated with this go term
        gene_set_in_target = gene_set[np.in1d(gene_set,target_genes)]
        b = gene_set_in_target.size        
        if b != 0:
            # calculate the enrichment probability as the cumulative sum of the tail end of a hypergeometric distribution
            # with parameters (N,B,n,b)
            n = target_genes.size
            num_iter = min(n,B)
            rng = np.arange(b,num_iter+1)
            probs.append(sum([np.exp(log_binomial(n,i)+log_binomial(N-n,B-i) - log_binomial(N,B)) for i in rng]))
        else:
            probs.append(1.0)
        
        #append associated genes to a list
        probs_genes.append(gene_set_in_target)
        
    probs = np.array(probs)    
    probs_genes = np.array(probs_genes)
    
    # adjust p value to correct for multiple testing
    fdr_q_probs = probs.size*probs / rankdata(probs,method='ordinal')
    
    # filter out go terms based on the FDR q value and p value thresholds
    filt = np.logical_and(fdr_q_probs<fdr_thresh,probs<p_thresh)
    enriched_goterms = goterms[filt]
    p_values = probs[filt]
    fdr_q_probs = fdr_q_probs[filt]    
    probs_genes=probs_genes[filt]
    
    # construct the Pandas DataFrame
    gns = []
    for i in probs_genes:
        gns.append(';'.join(i))
    gns = np.array(gns)
    enriched_goterms = pd.DataFrame(data=fdr_q_probs,index=enriched_goterms,columns=['fdr_q_value'])
    enriched_goterms['p_value'] = p_values
    enriched_goterms['genes'] = gns
    
    # sort in ascending order by the p value
    enriched_goterms = enriched_goterms.sort_values('p_value')   
    return enriched_goterms
