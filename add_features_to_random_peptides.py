import pandas as pd
input = '6000_random_spliced_9mers_with_cleavage_sites.csv'

def add_features_df(input):
    df = pd.read_csv(input,index_col=0)
    df = df[~df['sequence'].str.contains('X')]
    df['sequence1'] = [x[0:y] for x,y in zip(df['sequence'],df['length1'])]
    df['sequence2'] = [x[y:z] for x, y, z in zip(df['sequence'], df['length1'], df['length'])]
    df['n_term'] = [x[0] for x in df['sequence']]
    df['c_term'] = [x[-1] for x in df['sequence']]
    df['bind1'] = [x[y-1] for x,y in zip(df['sequence'],df['length1'])]
    df['bind2'] = [x[y] for x, y in zip(df['sequence'], df['length1'])]
    sequence_df = pd.DataFrame([list(x) for x in df['sequence']])
    sequence_df.columns = ['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9']
    
    cleavage_df = pd.DataFrame()
    for cleavage in ['cleavage1','cleavage2','cleavage3','cleavage4']:
        #cleavage_list = [list(x) for x in df[cleavage]]
        for i in range(0,5):
            cleavage_df[cleavage + "_" + str(i)] = pd.DataFrame([list(x) for x in df[cleavage]]).iloc[:,i]
    final_df = pd.concat([sequence_df,df[['length1','length2','length','n_term','c_term',
                                          'bind1','bind2','distance','reversed']],cleavage_df], axis=1)
    
    #df[['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9','seq_10','seq_11','seq_12']] = [list(x) for x in df['sequence']]
    final_df['output'] = 0
    final_df.to_csv('training_data_3_negative.csv')
    #pd.get_dummies(final_df.applymap(str),dummy_na=True).to_csv('nn_positive_input.csv')
    
add_features_df(input)