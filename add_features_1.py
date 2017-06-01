import csv
import pandas as pd
input = 'csv_files_parsed_peptides/all_spliced_peptides.csv'
output = 'training_data_3_positive.csv'
amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

def add_features(input,output):
    '''
    
    :param input: [sequence,protein,startpos1,endpos1,startpos2,endpos2]
    :param output: [sequence,protein,startpos1,endpos1,startpos2,endpos2,length1,length2,length,start1,end1,start2,end2]
    :return: 
    '''
    with open(input, 'r') as f, open(output,'w') as g:
        reader = csv.reader(f)
        writer = csv.writer(g)
        for row in reader:
            row[2:6] = [int(x) for x in row[2:6]]
            length1 = row[3] + 1 - row[2]
            length2 = row[5] + 1 - row[4]
            total_length = length1 + length2
            sequence = row[0]
            sequence1 = sequence[0:length1]
            sequence2 = sequence[length1:total_length]
            n_term = sequence1[0]
            bind1 = sequence1[length1-1]
            bind2 = sequence2[0]
            c_term = sequence[total_length-1]
            row.extend([length1,length2,total_length,n_term,bind1,bind2,c_term])
            writer.writerow(row)

            
def add_features_df(input,output):
    df_proteome = pd.read_csv('netchopoutput.csv')
    df_proteome.columns = ['position','aa','cleavage','probability','protein']
    df_proteome = df_proteome[['position','aa','protein']]
    df_proteome.set_index(['protein','position'])
    df = pd.read_csv(input)
    df.columns = ["sequence","protein","startpos1","endpos1","startpos2","endpos2"]
    df['distance'] = [min(abs(s1-e0),abs(s0-e1)) for s0,e0,s1,e1 in zip(df['startpos1'],df['endpos1'],df['startpos2'],df['endpos2'])]
    df['reversed'] = [s0>s1 for s0,s1 in zip(df['startpos1'],df['startpos2'])]
    df['length1'] = [x - y + 1 for x, y in zip(df['endpos1'],df['startpos1'])]
    df['length2'] = [x - y + 1 for x, y in zip(df['endpos2'], df['startpos2'])]
    df['length'] = [x + y for x, y in zip(df['length1'], df['length2'])]
    df['sequence1'] = [x[0:y] for x,y in zip(df['sequence'],df['length1'])]
    df['sequence2'] = [x[y:z] for x, y, z in zip(df['sequence'], df['length1'], df['length'])]
    df['n_term'] = [x[0] for x in df['sequence']]
    df['c_term'] = [x[-1] for x in df['sequence']]
    df['bind1'] = [x[y-1] for x,y in zip(df['sequence'],df['length1'])]
    df['bind2'] = [x[y] for x, y in zip(df['sequence'], df['length1'])]
    
    df_cleavages = pd.DataFrame()
    for site, position in zip(['cleavage1','cleavage3'],['startpos1','startpos2']):
        for i in range(5): # End residue is immediately after central resiude
            df_new = pd.DataFrame([[x-3+i,"sp|"+y+"|"] for x,y in zip(df[position],df['protein'])])
            df_new.columns = ['position','protein']
            col_name = site + "_" +str(i)
            df_cleavages[col_name] = pd.merge(df_new,df_proteome,how="left",on=['protein','position'])['aa']
    for site, position in zip(['cleavage2', 'cleavage4'], ['endpos1', 'endpos2']):
        for i in range(5): # This time the end residue is the central residue
            df_new = pd.DataFrame([[x-2+i,"sp|"+y+"|"] for x,y in zip(df[position],df['protein'])])
            df_new.columns = ['position','protein']
            col_name = site + "_" +str(i)
            df_cleavages[col_name] = pd.merge(df_new,df_proteome,how="left",on=['protein','position'])['aa']
    '''for index, row in df.iterrows():
        protein = "sp|"+ row['protein'] + "|"
        df_protein = df_proteome[(df_proteome['protein']==protein)]
        cleavage1 = df_protein[(df_protein['position']<=row['startpos1'])&(df_protein['position']>=row['startpos1']-2)].iloc[:,1].tolist()
        print(cleavage1)'''
    
    sequence_df = pd.DataFrame([list(x) for x in df['sequence']])
    sequence_df.columns = ['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9','seq_10','seq_11','seq_12']
    final_df = pd.concat([sequence_df,df[['length1','length2','length','n_term','c_term','bind1','bind2','distance','reversed']],df_cleavages], axis=1)
    #df[['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9','seq_10','seq_11','seq_12']] = [list(x) for x in df['sequence']]
    final_df['output'] = 1
    #pd.get_dummies(final_df.applymap(str),dummy_na=True).to_csv('nn_positive_input.csv')
    final_df.to_csv(output)
    
def add_features_df2(input):
    df = pd.read_csv(input)
    df.columns = ["sequence","protein","startpos1","endpos1","startpos2","endpos2"]
    df['distance'] = [min(abs(s1-e0),abs(s0-e1)) for s0,e0,s1,e1 in zip(df['startpos1'],df['endpos1'],df['startpos2'],df['endpos2'])]
    df['reversed'] = [s0>e1 for s0,e1 in zip(df['startpos1'],df['endpos2'])]
    df['length1'] = [x - y + 1 for x, y in zip(df['endpos1'],df['startpos1'])]
    df['length2'] = [x - y + 1 for x, y in zip(df['endpos2'], df['startpos2'])]
    df['length'] = [x + y for x, y in zip(df['length1'], df['length2'])]
    df['sequence1'] = [x[0:y] for x,y in zip(df['sequence'],df['length1'])]
    df['sequence2'] = [x[y:z] for x, y, z in zip(df['sequence'], df['length1'], df['length'])]
    df['n_term'] = [x[0] for x in df['sequence']]
    df['c_term'] = [x[-1] for x in df['sequence']]
    df['bind1'] = [x[y-1] for x,y in zip(df['sequence'],df['length1'])]
    df['bind2'] = [x[y] for x, y in zip(df['sequence'], df['length1'])]
    sequence_df = pd.DataFrame([list(x) for x in df['sequence']])
    sequence_df.columns = ['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9','seq_10','seq_11','seq_12']
    final_df = pd.concat([sequence_df,df[['length1','length2','length','n_term','c_term','bind1','bind2','distance','reversed']]], axis=1)
    #df[['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9','seq_10','seq_11','seq_12']] = [list(x) for x in df['sequence']]
    print(df['length'].value_counts())
    pd.get_dummies(final_df.applymap(str),dummy_na=True).to_csv('nn_positive_input.csv')
    
add_features_df(input,output)