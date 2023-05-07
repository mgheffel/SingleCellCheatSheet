#obs rewrite based on other columns 
adata.obs['temp_clust']=adata2.obs['leiden2'].apply(lambda i: clust_map[i] if i in clust_map else np.nan)
adata.obs['nL2']=np.where(adata.obs['temp_clust'].isnull(),adata.obs['nL2'],adata.obs['temp_clust'])

#Add obsm to obs
adata.obs[['umap_0','umap_1']]=adata.obsm['X_umap']

#harmony integration
sc.external.pp.harmony_integrate(adata, 'batch',max_iter_harmony=20)
sc.pp.neighbors(adata,n_pcs=20,n_neighbors=20,use_rep='X_pca_harmony')



      

