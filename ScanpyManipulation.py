#print all column options to rewrite
for i in sorted(adata.obs['leiden'].unique()):
    print("clust_map['"+i+"']='"+i+"'")
#obs rewrite based on other columns 
adata.obs['temp_clust']=adata.obs['leiden'].apply(lambda i: clust_map[i] if i in clust_map else np.nan)
adata.obs['col_name']=np.where(adata.obs['temp_clust'].isnull(),adata.obs['col_name'],adata.obs['temp_clust'])

for i in sorted(adata2.obs['panno'].unique()):
    print("adata2.obs['panno']=adata2.obs['panno'].replace('"+i+"','')")

#normalize features by global methylation levels
def normalize_dataframe_zero_one(df):
    min_values = df.min()
    max_values = df.max()
    df = (df - min_values) / (max_values - min_values)
    return df
adata.obs['mCG']=adata.obs['mCG/CG']-adata.obs['mCCC/CCC']
mcgdf=mcgdf.sub(adata[:,[g for g in adata.var.index if '_CG' in g]].obs['mCG'],axis=0)
mcgdf = normalize_dataframe_zero_one(mcgdf)
from scanpy import AnnData
adata=AnnData(mcgdf,obs=adata.obs,var=adata.var)

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

#downsample adata to all unique values of obs_col having same number of cells
def downsample(ad,obs_col,min_cluster_size=10):
    cdf=pd.DataFrame(ad.obs[obs_col].value_counts())
    cdf=cdf[cdf[cdf.columns[0]]>=min_cluster_size]
    ad2=ad[ad.obs[cdf.columns[0]].isin(cdf.index)].copy()
    keep_count=ad2.obs[obs_col].value_counts().min()
    print(keep_count)
    keep=[]
    for ct in ad2.obs[obs_col].unique():
        keep=keep+list(ad[ad.obs[obs_col]==ct].obs.sample(keep_count,random_state=0).index)
    ad2=ad[ad.obs.index.isin(keep)]
    ad2.obs[obs_col]=ad2.obs[obs_col].astype(str)
    ad2.obs[obs_col]=ad2.obs[obs_col].astype('category')
    ad2.var.index=[str(c) for c in ad2.var.index]
    return ad2

      

