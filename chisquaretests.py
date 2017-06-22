import pandas as pd
from scipy.stats import chi2_contingency
aa_dict = {'S': 0.083207136799280254, 'V': 0.05969329362753354, 'A': 0.070178284858173401, 'C': 0.022970256771492441, 'Q': 0.047658794128419044, 'P': 0.063103624877869691, 'N': 0.035881471493242882, 'T': 0.053571477460770595, 'G': 0.065702594916323465, 'M': 0.021322187589502981, 'X': 3.3560167684260995e-06, 'L': 0.099593771833428171, 'W': 0.012186050151078161, 'E': 0.071012696606281039, 'F': 0.03651443391901419, 'Y': 0.026648097884764452, 'H': 0.026268603041240058, 'I': 0.043406014352977189, 'K': 0.057270072888268402, 'D': 0.047394198701097866, 'R': 0.056413582082473761}

df_chi = []
for i in ['left_left','left','central','right','right_right']:
    df1 = pd.DataFrame.from_csv('splicing_rules/cleavage_frequencies_'+i+'_residue.csv',index_col=None,header=0)
    df2 = pd.DataFrame.from_csv('splicing_rules/cleavage_frequencies_'+i+'_residue_nonspliced.csv',index_col=None,header=0)
    df = pd.merge(df1,df2, on='aa')
    print(df)
    p = chi2_contingency(df[['count_x','count_y']].as_matrix().T)[1]
    df_chi.append([i,p])
df_chi= pd.DataFrame(df_chi,columns=['position','p-value'])
print(df_chi)
#df_chi.to_csv('splicing_rules/cleavage_positions_compared_first_and_second_spliced_peptides.csv')
    
for i in ['left_left','left','central','right','right_right']:
    #df1 = pd.DataFrame.from_csv('splicing_rules/cleavage_frequencies_'+i+'_residue_first_peptide.csv',index_col=None,header=0)
    df1 = pd.DataFrame.from_csv('splicing_rules/cleavage_frequencies_'+i+'_residue_second_peptide.csv',index_col=None,header=0)
    df_aa = pd.DataFrame.from_dict(aa_dict, orient='index')
    df_aa.index.name = 'aa'
    df_aa.reset_index(inplace=True)
    df_aa.columns = ['aa', 'Proteome Frequency']
    df2 = pd.merge(df1[['aa','count']],df_aa,on='aa')
    f_obs = df2['count'].as_matrix()
    f_exp = df2['Proteome Frequency'].as_matrix()*f_obs.sum()
    print(chisquare(f_obs,f_exp))
    print(f_exp,f_obs)
    p = chisquare(f_obs,f_exp)[1]
    df_chi.append([i,p])
df_chi= pd.DataFrame(df_chi,columns=['position','p-value'])
#df_chi.to_csv('splicing_rules/cleavage_positions_compared_nonspliced_and_proteome.csv')