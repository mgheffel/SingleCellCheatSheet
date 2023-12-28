#Normalize to zero one 
def normalize_col_zero_one(df,col):
    min_value = df[col].min()
    max_value = df[col].max()
    df[col] = df[col] - min_value
    df[col] = df[col] / (max_value - min_value)
    
#effeciency function to normalize an entire dataframe
def normalize_dataframe_zero_one(df):
    min_values = df.min()
    max_values = df.max()
    df = (df - min_values) / (max_values - min_values)
    return df
    
    
chunks_dict={}
chunks_dict['CA']=5
chunks_dict['Astro']=4
chunks_dict['GC']=4
chunks_dict['OPC']=2
chunks_dict['ODC']=2
for i in pdata.obs['dptL3'].unique():
    name=i
    if i in chunks_dict.keys():
        n_chunks=chunks_dict[i]
    else:
        n_chunks=1
    chunks = np.array_split(pdata[pdata.obs['dptL3']==name].obs.sort_values('inverted_dpt'), n_chunks)
    c=0
    for i in chunks:
        pdata.obs.loc[i.index, 'pdpt'] = name+'_'+str(c)
        c+=1
