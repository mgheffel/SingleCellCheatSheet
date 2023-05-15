#obs rewrite based on other columns 
adata.obs['temp_clust']=adata2.obs['leiden2'].apply(lambda i: clust_map[i] if i in clust_map else np.nan)
adata.obs['nL2']=np.where(adata.obs['temp_clust'].isnull(),adata.obs['nL2'],adata.obs['temp_clust'])

def load_obsm_to_obs(ad,key):
    ad.obs[[key+'_0',key+'_1']]=ad.obsm[key]
def load_obsm_from_obs(ad,key,obsm_key='X_umap'):
    ad.obsm[obsm_key]=ad.obs[[key+'_0',key+'_1']].to_numpy()
load_obsm_to_obs(adata2,'X_umap')
adata.obs['3C_X_umap_0']=adata2.obs['X_umap_0']
adata.obs['3C_X_umap_1']=adata2.obs['X_umap_1']
load_obsm_from_obs(adata,'3C_X_umap','X_umap')

#harmony integration
sc.external.pp.harmony_integrate(adata, 'batch',max_iter_harmony=20)
sc.pp.neighbors(adata,n_pcs=20,n_neighbors=20,use_rep='X_pca_harmony')



      

