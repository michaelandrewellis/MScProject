import pandas as pd
import numpy as np

def concatenate(start,end):
    finaldf = pd.read_excel("Liepe.xlsx",sheetname=start,header=0)
    finaldf = finaldf.iloc[:, [1, 3]]
    for i in range(start+1,end):
        print(i)
        df = pd.read_excel("Liepe.xlsx",sheetname=i,header=None)
        df = df.iloc[:,[1,3]]
        df.columns = finaldf.columns
        finaldf = pd.concat([finaldf,df],ignore_index=True)
    finaldf = finaldf.dropna()
    return(finaldf)
# There are still NaNs in here
df = concatenate(421,487)
df.to_csv('Fibroblasts_unspliced_peptides_1D_HCD.csv',index=False)

# 0-77
# 77-223
# 223-238
# 238-288
# 288-298
# 298-321
# 321-345
# 345-394
# 394-421
# 421-487

#df = pd.read_excel("Liepe.xlsx",sheetname=237,header=None)
#print(df)