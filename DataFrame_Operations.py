#Normalize to zero one 
def normalize_col_zero_one(df,col):
    min_value = pdf[col].min()
    max_value = pdf[col].max()
    df[col] = df[col] - min_value
    df[col] = df[col] / (max_value - min_value)
