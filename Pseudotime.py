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

#test root node options
for i in range(20):
    print(i)
    adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden']  == '5')[i]
    sc.tl.dpt(adata)
    sc.pl.umap(adata,color='dpt_pseudotime',vmax=.5)
    
#view root cells
for i in range(1,300):
    print(i)
    cell=adata[adata.obs['l53']=='RG-1_15'][i].obs.index[0]
    adata.obs['this_cell']=np.nan
    adata.obs.loc[cell,'this_cell']='1'
    adata.obs['dotsize']=1
    adata.obs.loc[cell,'dotsize']='20'
    sc.pl.umap(adata, color=['this_cell'],size=adata.obs['dotsize'])

#set root node and run pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden']  == '5')[8]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata,color='dpt_pseudotime')
sc.pl.umap(adata,color='dpt_pseudotime')

#invert paths for bidirectional paga plotting
invert_types = ['RG-2', 'OPC', 'ODC', 'Astro']
adata.obs['inverted_dpt'] = adata.obs['dpt_pseudotime'].copy()
adata.obs.loc[adata.obs['L3'].isin(invert_types), 'inverted_dpt'] =  adata.obs['dpt_pseudotime']*-1

               
               
paths={'Astro':['RG','RG-2','Astro'],
      'Olig':['RG','RG-2','ODC','OPC'],
      'UL':['RG','UL'],
      'DL':['RG','DL'],
      'CA':['RG','CA',],
      'DG':['RG','DG',],
      'ENT':['RG','ENT']}
path_pos_genes={}
path_neg_genes={}
for path in paths:
    adata2=adata[adata.obs['gL2'].isin([path])]
    path_pos_genes[path]=get_pseudotime_corr_vars(adata2,n_genes=20,negative=False)
    path_neg_genes[path]=get_pseudotime_corr_vars(adata2,n_genes=20,negative=True)
      
plt.rcParams['figure.figsize']=(10,4)
sc.pl.paga_path(adata3,adata3.obs.groupby('V2').mean().sort_values('dpt_pseudotime').index,reg_corp, groups_key='V2',
               color_map='seismic',color_maps_annotations={'dpt_pseudotime': 'viridis'},n_avg=100,show_node_names=True,
      

               ytick_fontsize=8)

###
#different size heatmap rows
#https://stackoverflow.com/questions/48459801/matplotlib-seaborn-control-line-row-height-of-heatmap
import numpy as np;np.random.seed(1)
import matplotlib.pyplot as plt

# get some data
a = np.random.rayleigh(3,111)
h,_ = np.histogram(a)
data = np.r_[[h]*10].T+np.random.rand(10,10)*11

# produce scaling for data
y = np.cumsum(np.append([0],np.sum(data, axis=1)))
x = np.arange(data.shape[1]+1)
X,Y = np.meshgrid(x,y)
# plot heatmap
im = plt.pcolormesh(X,Y,data)

# set ticks
ticks = y[:-1] + np.diff(y)/2
plt.yticks(ticks, np.arange(len(ticks)))
plt.xticks(np.arange(data.shape[1])+0.5,np.arange(data.shape[1]))
# colorbar
plt.colorbar(im)

plt.show()
###


###
###different colored heatmap rows
#https://stackoverflow.com/questions/68837536/how-to-use-a-different-colormap-for-different-rows-of-a-heatmap
# Reds
data1 = data.copy()
data1.loc[7] = float('nan')
ax = sns.heatmap(data1, annot=True, cmap="Reds")

# Greens
data2 = data.copy()
data2.loc[:6] = float('nan')
sns.heatmap(data2, annot=True, cmap="Greens")

###


###
custom pseudotime path genes chunk ct plot
df=pdata.to_df()
df['pseudotime']=pdata.obs['dpt_pseudotime']
df['inverted_dpt']=pdata.obs['inverted_dpt']
df['dptL3']=pdata.obs['pdpt']
df=pd.DataFrame(df.groupby('dptL3').mean())
df=df.transpose()


order=pdata.obs.groupby('pdpt').mean().sort_values('inverted_dpt').index
order
order=['ODC_0', 'ODC_1', 'OPC_0', 'OPC_1', 'Astro_0', 'Astro_1', 'Astro_2',
       'Astro_3', 'RGCs_G_0', 'RG-2_0', 'RGCs_0', 'RGCs_N_0', 'HBA2_0',
       'INPs_0', 'CA_0',  'CA_1', 'CA_2', 'CA_3', 'CA_4',
        'GC_0','GC_1','GC_2', 'GC_3']
df=df
df=df[order]
order

import seaborn as sns
plt.rcParams['figure.figsize']=(5,2)
pdf=df

pdf=pdf.transpose()
# pdf['pseudotime']=pdata.obs['inverted_dpt']
pdf=pdf[['CRYAB','MBP','FBXO32','DOCK10','SOX10',
         'CSPG4','PCDH15','OLIG1','OLIG2','PDGFRA',
         'DGKG','ALDH1A1',
         'FAM107A','HOPX',
         'COL5A3','ATP13A4','LRIG1','TNC','NES','PTN',
         
         'PAX6','SOX2','VIM','LUZP2','SPON1','EGFR',
         'GLI3','HES1','SFRP1','LIPG','FABP7','HBA2','SATB2',
         'SOX4','NEUROD2','SYT4',
         'EOMES','SOX11','TBR1',
         
         
         
         'NRG1',
         'MYT1L','GALNT18',
        'CSMD1',
        'RBFOX1','NELL2',
         
         'RYR3','DOCK9','CHRM3','KIAA0319','MEIS2',
        'SSTR2',
         'PROX1','TLL1','SEMA5A','CFAP43','pseudotime']]##cell
for col in pdf.columns:
    if 'dpt' in col or 'pseudotime' in col:
        continue
    normalize_col_zero_one(pdf,col)
# plt.rcParams['figure.figsize']=(8,12)
plt.rcParams['figure.figsize']=(6,6)
pdf=pdf.transpose()
pdf['b1']=np.nan
pdf['b2']=np.nan
pdf['b3']=np.nan
pdf['b4']=np.nan
pdf=pdf[['ODC_0', 'ODC_1', 'OPC_0', 'OPC_1','b1' ,'Astro_0', 'Astro_1', 'Astro_2',
       'Astro_3', 'RGCs_G_0', 'b3','RG-2_0', 'RGCs_0','RGCs_N_0',
         'b4',  'HBA2_0',
       'INPs_0', 'CA_0',  'CA_1', 'CA_2', 'CA_3', 'CA_4','b2',
        'GC_0','GC_1','GC_2', 'GC_3']]
pdf=pdf.transpose()
sns.set(font_scale=0.5)
sns.heatmap(pdf.transpose(), annot=False,cmap='viridis')
plt.savefig('genes_dpt.svg')
plt.show()
###


###Downsample

adata2=adata.copy()
keep=[]
keep=keep+list(adata2[adata2.obs['2T_L3']=='RG-1'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='RG-2'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='RG-CA'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='RG-UL'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='RG-DG'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='Exc-UL'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='Exc-CA-1'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='Exc-CA-3'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='Exc-DL'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='Exc-DL-2'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='Exc-ENT'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='RG-DG'].obs.sample(250,random_state=0).index)
keep=keep+list(adata2[adata2.obs['2T_L3']=='Oligo'].obs.sample(250,random_state=0).index)

print(adata2.shape)
adata2=adata2[adata2.obs.index.isin(keep)]