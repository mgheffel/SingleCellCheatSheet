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
               color_map='seismic',color_maps_annotations={'dpt_pseudotime': 'viridis'},n_avg=100,show_node_names=True,ytick_fontsize=8)

###
#multimodal compare paga plots
#create multimod adata for plotting
jdata=adata3.copy().transpose().concatenate(gsdf.copy().transpose()).transpose()
jdata=jdata[jdata.obs.index.isin(adata3.obs.index)]
jdata.var.index=[i[:-2] for i in jdata.var.index]
#optional check gene lengths
# df=pd.read_csv('../DATA/gencode.v33.gene.bed',sep='\t',header=None).set_index(3)
# df['len']=df[2]-df[1]
# df=df[~df.index.duplicated()]
# gsv=gsdf.var
# gsv=gsv[~gsv.index.duplicated()]
# gsv=gsv[gsv.index.isin(df.index)]
# gsv['len']=df['len']
# gsv=gsv.sort_values('len',ascending=False)
#plotting function
def plot_multimod(jdata,genes,order,obs_key,no_CH=False):
    rdf=[]
    jdf=jdata.to_df()
    jdf['b1']=.5
    count=0
    jgenes=[]
    for g in genes:
        if g in jdata.var.index and g+'_CG' in jdata.var.index :
            if no_CH or g+'_CH'in jdata.var.index:
                if gsv.loc[g][0]<100000:
                    try:
                        rdf.append((g,gsv.loc[g][0],'drop : len < 100,000'))
                    except:
                        rdf.append((g,np.nan,'drop : len < 100,000'))
                    print(g,' -drop : len < 100,000')
                    continue
            rdf.append((g,gsv.loc[g][0],''))
            jgenes.append(g)
            jgenes.append(g+'_CG')
            if not no_CH:
                jgenes.append(g+'_CH')
            jdf['b'+str(count)]=.5
            jgenes.append('b'+str(count))
            count+=1
        else:
            print(g,'else')
            try:
                rdf.append((g,gsv.loc[g][0],'drop : not in all mods'))
            except:
                rdf.append((g,np.nan,'drop : not in all mods'))
    jdata2=AnnData(jdf,obs=jdata.obs)
    print(jgenes)
    sc.pl.paga_path(jdata2,order,jgenes, groups_key=obs_key,
               color_map='seismic',color_maps_annotations={'dpt_pseudotime': 'viridis'},n_avg=100,
    normalize_to_zero_one=True,show_node_names=True,ytick_fontsize=6)
    return pd.DataFrame(rdf,columns=['gene','len','drop']).set_index('gene')
#call
plot_multimod(jdata,all_mod_genes,order=pathA,obs_key='panno',no_CH=True)
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
for ct in adata2.obs['sL3'].unique():
    keep=keep+list(adata2[adata2.obs['sL3']==ct].obs.sample(222,random_state=0).index)
print(adata2.shape)
data2=adata2[adata2.obs.index.isin(keep)]
print(adata2.shape)
