import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
sc.set_figure_params(figsize=(4,4))
matplotlib.rcParams['pdf.fonttype']=42
plt.rcParams['font.family']='Verdana'

#UMAP / draw_graph / scatter plot coordinate manipulation framework
adata.obsm['X_draw_graph_fa_mC']=adata.obsm['X_draw_graph_fa'][:, [1, 0]]
adata.obsm['X_pca_3C']=adata2.obsm['X_pca']
array=adata2.obsm['X_draw_graph_fa'][:, [1, 0]].copy()
array[:, 0] = -array[:, 0]
adata.obsm['X_draw_graph_fa_3C'] = array
array=adata2.obsm['X_umap'][:, [1, 0]].copy()
array[:, 0] = -array[:, 0]
adata.obsm['X_umap_3C']=array
#[:, [1, 0]]*-1
adata.obs['3C_dpt']=adata2.obs['dpt_pseudotime']


#expects 2d dataframe of fractions
from matplotlib.collections import PatchCollection
def dotplot_df(df,cmap='Reds'):
    xlabels=df.columns
    ylabels=df.head(40).index

    N = len(ylabels)
    M = len(xlabels)


    x, y = np.meshgrid(np.arange(M), np.arange(N))
    dfmax=df.max().max()
    df=np.clip(df,dfmax*-1,dfmax)
    s = df.head(40).to_numpy()
    c = df.head(40).to_numpy()

    fig, ax = plt.subplots()
    s=s*-1
    s=s-s.min()
    R = s/s.max()/2.5
    circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    col = PatchCollection(circles, array=c.flatten(), cmap=cmap)
    ax.add_collection(col)

    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.grid(which='minor')
    plt.xticks(rotation = 90)
    fig.colorbar(col)
    plt.show()
    return col
# dotplot_df(df.head(40))


#cell type donor ratio / number of cells 
ag='2T'
reg='HPC'

df2T=df[df['age_groups']==ag]
df2T=df2T[df2T['region']==reg]
df2T=df2T[~df2T['Sample'].isin(['20210316_UMB4267','20210316_UA1822'])]
df2T=pd.crosstab(df2T['L3'],df2T['Sample'])
sums=df2T.transpose().sum()
df2T['sum']=sums
df2T=df2T[df2T['sum']>=10]
sums=df2T.transpose().sum()
df2T=df2T[df2T.columns[:-1]]
df2T=df2T.transpose()
for c in df2T.columns:
    df2T[c]=df2T[c]/sums[c]
df2T=df2T.transpose()
plt.rcParams['figure.figsize']=(2,len(df2T)/4)
ax=df2T.plot(kind='barh',stacked=True,width=.9,color=sc.pl.palettes.vega_10[:len(df2T.columns)])
for c in ax.containers:
    labels = [str((v.get_width()*100).round(1))+'%' if v.get_width() > 0 else '' for v in c]
ax.legend(loc='center left',bbox_to_anchor=(1,.5))
plt.savefig(ag+'_'+reg+'_celltype_sample_comp.svg', bbox_inches='tight')
plt.show()

sdf=pd.DataFrame(sums)
sdf=np.log10(sdf)
ax=sdf.plot(kind='barh',width=.9,color=sc.pl.palettes.vega_10[:len(df2T)])
ax.set_xlabel('log10 cell count')
ax.axvline(x=2, linestyle='--', color='gray')
ax.legend(loc='center left',bbox_to_anchor=(1,.5))
plt.savefig(ag+'_'+reg+'_celltype_count.svg', bbox_inches='tight')
plt.show()
