#Normalize to zero one 
min_value = df[col].min()
max_value = df[col].max()
df[col] = df[col] - min_value
df[col] = df[col] / (max_value - min_value)

for col in ['mCH/CH','mCG/CG']:
    min_value = df[col].min()
    max_value = df[col].max()
    df[col] = df[col] - min_value
    df[col] = df[col] / (max_value - min_value)
