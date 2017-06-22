import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
aa_dict = {'S': 0.083207136799280254, 'V': 0.05969329362753354, 'A': 0.070178284858173401, 'C': 0.022970256771492441, 'Q': 0.047658794128419044, 'P': 0.063103624877869691, 'N': 0.035881471493242882, 'T': 0.053571477460770595, 'G': 0.065702594916323465, 'M': 0.021322187589502981, 'X': 3.3560167684260995e-06, 'L': 0.099593771833428171, 'W': 0.012186050151078161, 'E': 0.071012696606281039, 'F': 0.03651443391901419, 'Y': 0.026648097884764452, 'H': 0.026268603041240058, 'I': 0.043406014352977189, 'K': 0.057270072888268402, 'D': 0.047394198701097866, 'R': 0.056413582082473761}
amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

spliced_peptides = 'training_data_3_positive.csv'
nonspliced_peptides = 'csv_files_parsed_peptides/nonspliced_peptides_with_cleavage_sites.csv'

df = pd.read_csv(spliced_peptides, index_col=0)
df = df[df['length']==9].drop(['seq_10','seq_11','seq_12'],axis=1)
df = df.dropna()
df.drop_duplicates(inplace=True)

def get_distances(df):
    df = df[['distance','reversed']]
    df = df[df['reversed']==False]
    df['distance'].value_counts().sort_index().plot(kind='bar')
    plt.show()

def get_pairs(df):
    df=pd.DataFrame.from_csv('csv_files_parsed_peptides/all_unspliced_peptides.csv',index_col=None,header=None)
    df.columns = ['sequence','protein','start','end','length']
    pairlist = []
    for i in df['sequence']:
        for j in range(8):
            pair=i[j:j+2]
            pairlist.append(pair)
    df_pair=pd.Series(pairlist)
    df_pair = df_pair.value_counts()
    df_pair.to_frame().to_csv('splicing_rules/pair_count_nonspliced.csv')
    pair_table = []
    '''
    # limit to peptides not affected by achor preferences
    df_binding = df
    total = df_binding.shape[0]
    for i in amino_acids:
        for j in amino_acids:
            pair_prob = aa_dict[i]*aa_dict[j]
            try:
                pair_count = df_binding.groupby(['bind1','bind2']).count().loc[i,j][0]
            except:
                pair_count=0
            pair_freq = pair_count/total
            pair_table.append([i+j,pair_count,pair_freq,pair_prob])
            
    df = pd.DataFrame(pair_table)
    df.columns = ['pair','count','observed','expected']
    df[['pair','count']].to_csv('binding_pair_count_all_6_6_2017.csv')'''
    #df['p-value'] = [stats.binom_test(x,total,y,alternative='greater') for x,y, in zip(df['count'],df['expected'])]
   # df.sort_values(by='p-value').to_csv('anchor_points_2_9_excl_binding.csv',header=True,index=False)

def get_residues(df):
    aa_table = []
    df_cleavage = df
    total = df_cleavage.shape[0]
    for i in amino_acids:
        aa_prob = aa_dict[i]
        try:
            aa_count = df_cleavage.groupby(['bind2']).count().loc[i][0]
        except:
            aa_count=0
        aa_freq = aa_count/total
        aa_table.append([i,aa_count,aa_freq,aa_prob])
    df = pd.DataFrame(aa_table)
    df.columns = ['aa','count','observed','expected']
    df['p-value'] = [stats.binom_test(x,total,y,alternative='greater') for x,y, in zip(df['count'],df['expected'])]
    df.sort_values(by='p-value').to_csv('splicing_rules/bind_frequencies_right_residue_latest.csv',header=True,index=False)
    

def get_unspliced_residues():
    df = pd.DataFrame.from_csv('csv_files_parsed_peptides/unspliced_sequences.csv',header=0,index_col=None)
    df_cleavage = df
    total = df_cleavage.shape[0]
    for j in df.columns:
        aa_table = []
        for i in amino_acids:
            aa_prob = aa_dict[i]
            try:
                aa_count = df_cleavage.groupby([j]).count().loc[i][0]
            except:
                aa_count=0
            aa_freq = aa_count/total
            aa_table.append([i,aa_count,aa_freq,aa_prob])
        df = pd.DataFrame(aa_table)
        df.columns = ['aa','count','observed','expected']
        df['p-value'] = [stats.binom_test(x,total,y,alternative='greater') for x,y, in zip(df['count'],df['expected'])]
        df.sort_values(by='p-value').to_csv('splicing_rules/residue_frequency_'+j+'_unspliced.csv',header=True,index=False)

def find_anchor_residues(df):
    pair_table = []
    df_binding = df
    total = df_binding.shape[0]
    df_spliced_freq = pd.DataFrame()
    for pos in range(1,10):
        for i in amino_acids:
            count = df_binding.groupby('seq_'+str(pos)).count().loc[i][0]
            expected_freq = count/total
            df_spliced_freq = df_spliced_freq.append([[pos,i,expected_freq]],ignore_index=True)
    df_spliced_freq.columns = ['position','aa','frequency']
    triplet_list = []
    for index,row in df_binding.iterrows():
        for pos1 in range(1, 8):
            for pos2 in range(pos1 + 1, 9):
                for pos3 in range(pos2 + 1, 10):
                    triplet_list.append([pos1,pos2,pos3,row['seq_'+str(pos1)],row['seq_'+str(pos2)],row['seq_'+str(pos3)]])
    df = pd.DataFrame(triplet_list)
    df.columns = ['pos1','pos2','pos3','residue1','residue2','residue3']
    df = df.groupby(['pos1','pos2','pos3','residue1','residue2','residue3']).size().reset_index().rename(columns={0:'count'})
    df = pd.merge(df,df_spliced_freq,left_on=['pos1','residue1'],right_on=['position','aa'])
    df = pd.merge(df, df_spliced_freq, left_on=['pos2', 'residue2'], right_on=['position', 'aa'])
    df = pd.merge(df, df_spliced_freq, left_on=['pos3', 'residue3'], right_on=['position', 'aa'])
    df = df[['pos1','pos2','pos3','residue1','residue2','residue3','count','frequency_x','frequency_y','frequency']]
    df.columns = ['pos1','pos2','pos3','residue1','residue2','residue3','count','prob1','prob2','prob3']
    df['expected']= df['prob1'] * df['prob2'] * df['prob3']
    df.to_csv('all_anchor_triples_test.csv', header=True, index=False)
    df['p-value'] = [stats.binom_test(x, total, y, alternative='greater') for x, y, in zip(df['count'], df['expected'])]
    df.sort_values(by='p-value').to_csv('all_anchor_triples_using_actual_frequencies.csv', header=True, index=False)
    
    
def binding_overrep():
    df_spliced = pd.DataFrame.from_csv('splicing_rules/residue_frequencies_all_positions_spliced_and_unspliced.csv',index_col=None)
    df_spliced = df_spliced[df_spliced['type']=='spliced'].drop('type',axis=1)
    df_spliced = df_spliced[['aa','count']].groupby('aa').sum()
    df_spliced['expected_freq'] = [df_spliced.loc[x][0]/df_spliced.sum()[0] for x in df_spliced.index]
    print(df_spliced)
    df_spliced.reset_index(inplace=True)
    df_binding_left = pd.DataFrame.from_csv('splicing_rules/bind_frequencies_right_residue.csv',index_col=None)[['aa','count']]
    df = pd.merge(df_spliced[['aa','expected_freq']],df_binding_left, on='aa')
    total = df['count'].sum()
    df['observed_freq'] = df['count']/total
    df['p-value'] = [stats.binom_test(x, total, y) for x, y, in zip(df['count'], df['expected_freq'])]
    df.sort_values(by='p-value').to_csv('binding_frequency_right_relative_to_non_binding_in_spliced_peptides.csv', header=True, index=False)

def binding_pair_overrep():
    df_spliced = pd.DataFrame.from_csv('splicing_rules/residue_frequencies_all_positions_spliced_and_unspliced.csv',index_col=None)
    df_spliced = df_spliced[df_spliced['type']=='spliced'].drop('type',axis=1)
    df_spliced = df_spliced[['aa','count']].groupby('aa').sum()
    df_spliced['expected_freq'] = [df_spliced.loc[x][0]/df_spliced.sum()[0] for x in df_spliced.index]
    print(df_spliced)
    df_spliced.reset_index(inplace=True)
    df_pairs = pd.DataFrame.from_csv('splicing_rules/binding_pair_frequency_compared_to_proteome.csv', index_col=None)[['pair', 'count']]
    df_pairs['bind1'] = [list(str(x))[0] for x in df_pairs['pair']]
    df_pairs['bind2'] = [list(str(x))[1] for x in df_pairs['pair']]
    df_pairs.drop('pair',inplace=True,axis=1)
    df = pd.merge(df_spliced[['aa','expected_freq']],df_pairs, right_on='bind1', left_on='aa')
    df = pd.merge(df_spliced[['aa', 'expected_freq']], df, right_on='bind2', left_on='aa')
    print(df)
    total = df['count'].sum()
    df['observed_freq'] = df['count']/total
    df['p-value'] = [stats.binom_test(x, total, y*z) for x, y, z in zip(df['count'], df['expected_freq_x'],df['expected_freq_y'])]
    df.sort_values(by='p-value').to_csv('binding_pair_frequency_relative_to_non_binding_in_spliced_peptides.csv', header=True, index=False)
    print(df)
    
def proline_analysis(df):
    df = df[(df['seq_2']=='R')]
    df = df.groupby('length1').count()['output']
    df.plot(kind='bar',use_index=False)
    plt.show()

def get_gap(df):
    df = df[df['distance']==2]
    df = df[df['reversed']==False]
    df = df['cleavage3_2'].value_counts().to_frame()
    df['frequency'] = df['cleavage3_2']/df['cleavage3_2'].sum()
    df_aa = pd.DataFrame.from_dict(aa_dict,orient='index')
    df = pd.merge(df,df_aa,left_index=True,right_index=True)
    df.columns = ['f_obs','frequency','proteome_frequency']
    df['f_exp'] = df['proteome_frequency']*df['f_obs'].sum()
    df[['f_exp','f_obs']].plot(kind='bar')
    plt.show()
    
    
get_gap(df)
#get_unspliced_residues()
#binding_pair_overrep()
#find_anchor_residues(df)
#proline_analysis(df)