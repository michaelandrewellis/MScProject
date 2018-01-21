import pandas as pd
import numpy as np
allele_list = ['A0301', 'A0101', 'B0702', 'B2705', 'C0702', 'C0202']


for a in allele_list:
    file = 'all_invented_peptides_50000_'+a+'.out'
    with open(file, 'r') as f:
        text = f.read().splitlines()
        text = text[32:-1]
        table = [line.split() for line in text]
        new_df = pd.DataFrame(table)
        new_df.columns = [0,'allele','sequence','PEPLIST','affinity']
        allele=new_df.allele[1]
        new_df.set_index('sequence',inplace=True)
        new_df = new_df['affinity']
        new_df.name=allele
        try:
            df = pd.concat([df, new_df],axis=1)
        except:
            df = new_df
            
# Transform to nM            
df = np.exp(((1-df.astype(float))*np.log(50000)))

# Save binding predictions
df.to_csv('all_invented_peptides_50000_with_GRLCL_alleles.csv',index=True)   

# Save list of predicted binders
df_binders = df[(df<500).T.any()].index.to_series().reset_index(drop=True)
df_binders.rename('sequence',inplace=True)
df_binders.to_csv('all_predicted_binders_50000_GRLCL_invented.csv')

# Create training set of predicted binders
train = pd.DataFrame.from_csv('training_data_with_cleavage_negative_50000.csv')
train['sequence'] = train.iloc[:,-17:-8].sum(axis=1)
train.set_index('sequence',inplace=True)
binders = set(df_binders)
seqs_to_drop = []
for sequence,row in train.iterrows():
    if sequence not in binders:
        seqs_to_drop.append(sequence)
seqs_to_drop = list(set(seqs_to_drop))
train.drop(seqs_to_drop,inplace=True)
train.reset_index(drop=True,inplace=True)
train.to_csv('training_data_negative_binders_50000_GRLCL.csv')



