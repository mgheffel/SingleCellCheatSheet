#print all column options to rewrite
for i in sorted(adata.obs['leiden'].unique()):
    print("clust_map['"+i+"']='"+i+"'")
#obs rewrite based on other columns 
adata.obs['temp_clust']=adata.obs['leiden'].apply(lambda i: clust_map[i] if i in clust_map else np.nan)
adata.obs['col_name']=np.where(adata.obs['temp_clust'].isnull(),adata.obs['col_name'],adata.obs['temp_clust'])

for i in sorted(adata2.obs['panno'].unique()):
    print("adata2.obs['panno']=adata2.obs['panno'].replace('"+i+"','')")

def load_obsm_to_obs(ad,key):
    ad.obs[[key+'_0',key+'_1']]=ad.obsm[key]
def load_obsm_from_obs(ad,key,obsm_key='X_umap'):
    ad.obsm[obsm_key]=ad.obs[[key+'_0',key+'_1']].to_numpy()
load_obsm_to_obs(adata2,'X_umap')
adata.obs['3C_X_umap_0']=adata2.obs['X_umap_0']
adata.obs['3C_X_umap_1']=adata2.obs['X_umap_1']
load_obsm_from_obs(adata,'3C_X_umap','X_umap')

#invert umap along diagonal axis (swap X,Y)
adata.obs['mC_X_draw_graph_faF_0']=adata.obs['mC_X_draw_graph_fa_1']
adata.obs['mC_X_draw_graph_faF_1']=adata.obs['mC_X_draw_graph_fa_0']

#harmony integration
sc.external.pp.harmony_integrate(adata, 'batch',max_iter_harmony=20)
sc.pp.neighbors(adata,n_pcs=20,n_neighbors=20,use_rep='X_pca_harmony')

#remove 100kb bin variables
gvar=[]
for g in adata.var.index:
    try: int(g.split('_')[0])
    except: gvar.append(g)



      

