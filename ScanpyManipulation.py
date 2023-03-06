#obs rewrite based on other columns 
adata.obs['temp_clust']=adata2.obs['leiden2'].apply(lambda i: clust_map[i] if i in clust_map else np.nan)
adata.obs['nL2']=np.where(adata.obs['temp_clust'].isnull(),adata.obs['nL2'],adata.obs['temp_clust'])


#DPT Pseudotime

def get_differential_genes(ad,groupby='L3',num_genes=20):
    sc.tl.rank_genes_groups(ad,groupby=groupby)
    hv=[]
    for i in range(num_genes):
        for i in ad.uns['rank_genes_groups']['names'][i]:
            if i not in hv: hv.append(i)
        if len(hv)>=num_genes: break
                
    return list(set(hv))
def get_pseudotime_corr_vars(ad,n_genes=10,negative=False):
    xdf=ad.to_df()
    xdf['dpt']=ad.obs['dpt_pseudotime']
    corrs_list=[]
    for c in xdf.columns[:-1]:
        corrs_list.append(xdf[c].corr(xdf['dpt']))
    cdf=pd.DataFrame(corrs_list,index=xdf.columns[:-1])
    cdf=cdf.rename(columns={0:'corr'})
    if negative:
        return list(cdf.sort_values('corr').head(n_genes).index)
    else:
        return list(cdf.sort_values('corr',ascending=False).head(n_genes).index)

      
sc.tl.paga(adata2, groups='leiden')
sc.pl.paga(adata2,color=['lp25'],threshold=0)
sc.tl.draw_graph(adata2, init_pos='paga')
sc.pl.draw_graph(adata2,color=['lp25'])
