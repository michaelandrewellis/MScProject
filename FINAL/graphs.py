import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save
import pandas as pd
import numpy as np
import pyfaidx as pf
import itertools

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

aa_cats = {'S': 'N', 'V': 'H', 'A': 'H', 'C': 'H', 'Q': 'N', 'P': 'P', 'N': 'N', 'T': 'N', 'G': 'H', 'M': 'H', 'L': 'H', 'W': 'W', 'E': 'E', 'F': 'W', 'Y': 'W', 'H': 'K', 'I': 'H', 'K': 'K', 'D':'E', 'R': 'K'}
aa_dict = {'S': 0.083207136799280254, 'V': 0.05969329362753354, 'A': 0.070178284858173401, 'C': 0.022970256771492441, 'Q': 0.047658794128419044, 'P': 0.063103624877869691, 'N': 0.035881471493242882, 'T': 0.053571477460770595, 'G': 0.065702594916323465, 'M': 0.021322187589502981, 'L': 0.099593771833428171, 'W': 0.012186050151078161, 'E': 0.071012696606281039, 'F': 0.03651443391901419, 'Y': 0.026648097884764452, 'H': 0.026268603041240058, 'I': 0.043406014352977189, 'K': 0.057270072888268402, 'D': 0.047394198701097866, 'R': 0.056413582082473761}
aa_dict_bind1 = {'G':0.1048723897911833,'L':0.10409899458623356,'A':0.10177880897138437,'S':0.08553750966744006,'P':0.08043310131477185,'R':0.07037896365042537,'V':0.06867749419953596,'K':0.05645784996133024,'T':0.048414539829853054,'I':0.04501160092807425,'E':0.04501160092807425,'Q':0.032946635730858466,'D':0.031245166279969063,'F':0.027532869296210363,'N':0.026295436968290797,'H':0.018716163959783448,'Y':0.017633410672853827,'M':0.01577726218097448,'C':0.014849187935034803,'W':0.004331013147718484}
aa_dict_spliced = {'L':0.13983778250255927,'A':0.08452634065674462,'P':0.08169147176943066,'S':0.07595873690841799,'R':0.07502952988424286,'G':0.07487203716828096,'V':0.07057248602252146,'K':0.06789510985116938,'I':0.05630364595637452,'T':0.04704307425781558,'E':0.043688479407827385,'Q':0.033577447043074256,'D':0.029466887156469012,'F':0.028962910465390977,'N':0.025639814158595165,'Y':0.02009607055673675,'H':0.01737144657059611,'M':0.011733207339160563,'C':0.01094574375935113,'W':0.0047877785652413575}
aa_dict_spliced_no_bind1 = {'A':0.08257057688935648,'C':0.010503243906715764,'D':0.029265298965456776,'E':0.043538488514816766,'F':0.02912502191828862,'G':0.07147115553217605,'H':0.017219007539891286,'I':0.057583727862528494,'K':0.0691916535156935,'L':0.14388918113273716,'M':0.011274767666140629,'N':0.025565491846396633,'P':0.08183412239172365,'Q':0.03364895668946169,'R':0.07555672453094862,'S':0.07487287392600386,'T':0.046887603015956514,'V':0.07078730492723129,'W':0.00483955812730142,'Y':0.02037524110117482}
aa_dict_hydrophobicity = {'S': 'ζ', 'V': 'Φ', 'A': 'Φ', 'C': 'Φ', 'Q': 'ζ', 'P': 'Φ', 'N': 'ζ', 'T': 'ζ', 'G': 'Φ', 
                          'M': 'Φ', 'L': 'Φ', 'W': 'Φ', 'E': 'ζ', 'F': 'Φ', 'Y': 'Φ', 'H': 'ζ', 'I': 'Φ', 'K': 'ζ', 
                          'D':'ζ', 'R': 'ζ'}
aa_dict_charge = {'S': 'P', 'V': 'N', 'A': 'N', 'C': 'P', 'Q': 'P', 'P': 'N', 'N': 'P', 'T': 'P', 'G': 'N', 
                          'M': 'N', 'L': 'N', 'W': 'N', 'E': '-', 'F': 'N', 'Y': 'P', 'H': '+', 'I': 'N', 'K': '+', 
                          'D':'-', 'R': '+'}
possible_pairs = {i[0]+i[1] for i in set(itertools.product(aa_dict.keys(),repeat=2))}
possible_triples = {i[0]+i[1]+i[2] for i in set(itertools.product(aa_dict.keys(),repeat=3))}

def plot_9mer_sequence(cell_line):
    df = pd.DataFrame.from_csv('all_peptides.csv')
    df = df[df['peptide'].str.len() == 9]
    df = df[df['cell_line']==cell_line]
    df.reset_index(inplace=True)
    sequence_df = pd.DataFrame([list(x) for x in df['peptide']])
    df = pd.concat([sequence_df,df],axis=1)
    spliced_total = df.groupby('type').count().iloc[0,0]
    unspliced_total = df.groupby('type').count().iloc[1,0]
    f, axes = plt.subplots(3, 3, sharey=True)
    for i in range(3):
        for j in range(3):
            axe = axes[i, j]
            k = 3 * i + j + 1
            spliced_series = df[df['type'] == 'spliced'].groupby(k - 1)[k - 1].count().apply(
                lambda x: x / spliced_total)
            unspliced_series = df[df['type'] == 'unspliced'].groupby(k-1)[k-1].count().apply(lambda x: x/unspliced_total)
            df_plot = pd.concat([spliced_series,unspliced_series],axis=1)
            ax = df_plot.plot(ax=axe, kind='bar')
            #ax.set_xticks(df_plot.index)
    plt.show()

def plot_9mer_sequence2(cell_line=None,position=1):
    df = pd.DataFrame.from_csv('all_peptides.csv')
    df = df[df['peptide'].str.len() == 9]
    df.drop_duplicates(subset=['peptide'],inplace=True)
    if cell_line is not None:
        df = df[df['cell_line']==cell_line]
    df.reset_index(inplace=True)
    df['aa'] = df.peptide.str[int(position)-1]
    spliced = df[df.type=='spliced']['aa'].value_counts(normalize=True).rename('spliced')
    unspliced = df[df.type=='unspliced']['aa'].value_counts(normalize=True).rename('unspliced')
    pd.concat([spliced,unspliced],axis=1).sort_values('unspliced').plot(kind='bar',color=['y','b'],edgecolor='k')
    plt.xlabel('Residue')
    plt.ylabel('Frequency')
    plt.yticks(np.arange(0,0.41,0.1))
    plt.xticks(rotation='horizontal')
    plt.title('Position 1')
    tikz_save('writeup/images/positions/9mers_'+str(cell_line)+'_'+str(position)+'.tex',
              figureheight = '\\figureheight',
              figurewidth = '\\figurewidth')

def plot_all_9mers(cell_line=None):
    for i in range(1,10):
        plot_9mer_sequence2(cell_line=cell_line,position=i)
    
def generate_binding_pairs(pair_type='centre'):
    # Compare binding pairs to expected pairs from unspliced peptides excluding anchor residues
    df_spliced = pd.DataFrame.from_csv('all_spliced_peptides_with_cell_lines.csv')
    df_spliced.drop('cell_line',inplace=True,axis=1)
    df_spliced.drop_duplicates(subset=['peptide'],keep=False,inplace=True)
    df_spliced['length1'] = df_spliced.end1-df_spliced.start1+1
    df_spliced['length2'] = df_spliced.end2-df_spliced.start2+1
    df_unspliced = pd.DataFrame.from_csv('all_unspliced_peptides_with_cell_lines.csv')
    pairs = []
    for peptide in df_unspliced.peptide:
        non_anchor = peptide[2:-1]
        pairs.extend([non_anchor[i:i+2] for i in range(0,5)]) #len(non_anchor)-2
        df_pairs_prob = pd.Series(pairs).value_counts()
    if pair_type == 'centre':
        df_spliced = df_spliced[(df_spliced.length1>2)&(df_spliced.length2>1)]
        df_binding_pairs = df_spliced.apply(lambda x: x[0][x[6]-1:x[6]+1],axis=1).value_counts()
    if pair_type == 'right':
        df_spliced = df_spliced[(df_spliced.length1>1)&(df_spliced.length2>2)]
        df_binding_pairs = df_spliced.apply(lambda x: x[0][x[6]:x[6]+2],axis=1).value_counts()
    if pair_type == 'left':
        df_spliced = df_spliced[(df_spliced.length1 > 3) & (df_spliced.length2 > 0)]
        df_binding_pairs = df_spliced.apply(lambda x: x[0][x[6]-2:x[6]], axis=1).value_counts()  
    df_binding_pairs = pd.concat([df_binding_pairs,df_pairs_prob],axis=1)
    df_binding_pairs.columns = ['count','unspliced_frequency']
    df_binding_pairs.to_csv('binding_pair_' + pair_type + '.csv')
    print(df_binding_pairs)
    
def plot_single_residue_between_peptides():
    fasta = pf.Fasta('../reviewedproteome.fasta',key_function = lambda x: x.split('|')[1])
    df = pd.DataFrame.from_csv('all_spliced_peptides_with_cell_lines.csv')
    df = df[((df['start2']-df['end1']) == 2)]
    nonambiguous = set(pd.DataFrame.from_csv('all_spliced_peptides_nonambiguous.csv')['0'])
    for item in df.index:
        if item not in nonambiguous:
            df = df.drop(item)
    df['single'] = [fasta[x][int(y)].seq for x, y in zip(df['protein'], df['start2'])]
    df['single'].value_counts(normalize=True).plot(kind='bar')
    plt.show()
    
def netchop_cleavage_predictions_spliced(type):
    # Doesn't show much. Perhaps only that cleavage site 4 looks most like a normal cleavage site. Maybe look at mean prediction scores.
    df = pd.DataFrame.from_csv('all_spliced_peptides_with_cell_lines.csv')
    # Remove duplicates and ambiguous peptides
    df.drop('cell_line',axis=1,inplace=True)
    df.drop_duplicates(subset=['peptide'],keep=False,inplace=True)
    nonambiguous = set(pd.DataFrame.from_csv('all_spliced_peptides_nonambiguous.csv')['0'])
    for item in df.peptide:
        if item not in nonambiguous:
            df = df[df.peptide!=item]
    df[['cleavage1','cleavage3']] = df[['start1','start2']].applymap(lambda x:str(int(x)-1))
    df[['cleavage2', 'cleavage4']] = df[['end1', 'end2']].applymap(lambda x:str(int(x)))
    df_chop = pd.DataFrame.from_csv('../netchopoutput20s.csv',index_col=None)
    df_chop.columns = ['position','aa','cleavage_orig','cleavage','protein'] # 'prob' changed to 'cleavage'
    df_chop.dropna(inplace = True)
    df_chop['protein'] = df_chop['protein'].apply(lambda x: x.split('|')[1])
    if type == 'predictions':
        series_cleavage1 = pd.merge(df, df_chop[['protein', 'position', 'cleavage_orig']], left_on=['protein', 'cleavage1'],
                                    right_on=['protein', 'position'], how='left')['cleavage_orig']
        series_cleavage2 = pd.merge(df, df_chop[['protein', 'position', 'cleavage_orig']], left_on=['protein', 'cleavage2'],
                                    right_on=['protein', 'position'], how='left')['cleavage_orig']
        series_cleavage3 = pd.merge(df, df_chop[['protein', 'position', 'cleavage_orig']], left_on=['protein', 'cleavage3'],
                                    right_on=['protein', 'position'], how='left')['cleavage_orig']
        series_cleavage4 = pd.merge(df, df_chop[['protein', 'position', 'cleavage_orig']], left_on=['protein', 'cleavage4'],
                                    right_on=['protein', 'position'], how='left')['cleavage_orig']
        pd.DataFrame([['cleavage1',series_cleavage1.value_counts(normalize=True)['S']],
                      ['cleavage2',series_cleavage2.value_counts(normalize=True)['S']],
                      ['cleavage3', series_cleavage3.value_counts(normalize=True)['S']],
                      ['cleavage4', series_cleavage4.value_counts(normalize=True)['S']]]).plot(kind='bar')
    elif type == 'scores':
        series_cleavage1 = pd.merge(df, df_chop[['protein', 'position', 'cleavage']], left_on=['protein', 'cleavage1'],
                                    right_on=['protein', 'position'], how='left')['cleavage']
        series_cleavage2 = pd.merge(df, df_chop[['protein', 'position', 'cleavage']], left_on=['protein', 'cleavage2'],
                                    right_on=['protein', 'position'], how='left')['cleavage']
        series_cleavage3 = pd.merge(df, df_chop[['protein', 'position', 'cleavage']], left_on=['protein', 'cleavage3'],
                                    right_on=['protein', 'position'], how='left')['cleavage']
        series_cleavage4 = pd.merge(df, df_chop[['protein', 'position', 'cleavage']], left_on=['protein', 'cleavage4'],
                                    right_on=['protein', 'position'], how='left')['cleavage']
        pd.DataFrame([['cleavage1', series_cleavage1.mean()],
                      ['cleavage2', series_cleavage2.mean()],
                      ['cleavage3', series_cleavage3.mean()],
                      ['cleavage4', series_cleavage4.mean()]]).plot(kind='bar')
    plt.show()
    
def netchop_cleavage_predictions_unspliced():
    # Doesn't show much. Perhaps only that cleavage site 4 looks most like a normal cleavage site. Maybe look at mean prediction scores.
    df = pd.DataFrame.from_csv('all_unspliced_peptides_with_cell_lines.csv')
    df.columns = ['sequence','protein','start1','end1','cell_line']
    df[['cleavage1']] = df[['start1']].applymap(lambda x:str(int(x)-1))
    df[['cleavage4']] = df[['end1']].applymap(lambda x:str(int(x)))
    df_chop = pd.DataFrame.from_csv('../netchopoutput20s.csv',index_col=None)
    df_chop.columns = ['position','aa','cleavage_orig','cleavage','protein'] # 'prob' changed to 'cleavage'
    df_chop.dropna(inplace = True)
    df_chop['protein'] = df_chop['protein'].apply(lambda x: x.split('|')[1])
    series_cleavage1 = pd.merge(df,df_chop[['protein','position','cleavage_orig']],left_on=['protein','cleavage1'],right_on=['protein','position'],how='left')['cleavage_orig']
    series_cleavage4= pd.merge(df,df_chop[['protein','position','cleavage_orig']],left_on=['protein','cleavage4'],right_on=['protein','position'],how='left')['cleavage_orig']
    pd.DataFrame([['cleavage1',series_cleavage1.value_counts(normalize=True)['S']],
                  ['cleavage4', series_cleavage4.value_counts(normalize=True)['S']]]).plot(kind='bar')
    '''pd.DataFrame([['cleavage1', series_cleavage1.mean()],
                  ['cleavage4', series_cleavage4.mean()]]).plot(kind='bar')'''
    plt.show()
    
def peptide_to_cleavage_comparison(cell_line=''):
    '''
    Compare cleaved residues to residues in peptide splicing partner. This is to investigate the possibility of splicing
    occuring via tranpeptidation?
    '''
    proteins = pf.Fasta('../reviewedproteome.fasta', key_function=lambda x: x.split('|')[1])
    df = pd.DataFrame.from_csv('all_spliced_peptides_with_cell_lines.csv')
    if cell_line:
        df = df[df['cell_line']==cell_line]
    df.drop('cell_line',axis=1,inplace=True)
    df.drop_duplicates(subset=['peptide'],keep=False,inplace=True)
    nonambiguous = set(pd.DataFrame.from_csv('all_spliced_peptides_nonambiguous.csv')['0'])
    for item in df.peptide:
        if item not in nonambiguous:
            df = df[df.peptide!=item]
    match_counts = [0]*8
    total_counts = [0]*8
    pairs = []
    cleavage_bind1_bind2 = []
    bind_right_triple=[]
    for row in df.iterrows():
        protein = row[1]['protein']
        end1 = row[1]['end1']
        end2 = row[1]['end2']
        start1 = row[1]['start1']
        start2 = row[1]['start2']
        length1 = end1-start1+1
        length2 = end2 - start2 +1
        distance = max(start1-end2,start2-end1)
        try:
            protein_length = proteins[protein].__len__()
        except:
            continue
        try:
            cleavage2 = proteins[protein][end1:min(end1+4,proteins[protein].__len__())].seq
            cleavage3 = proteins[protein][max(0,start2-5):start2-1].seq
        except KeyError:
            print('KeyError')
            continue
        except ValueError:
            print('protein: ',protein,proteins[protein].__len__())
            continue
        binding1 = row[1]['peptide'][max(length1-4,0):length1]
        binding2 = row[1]['peptide'][length1:min(length1+length2,length1+4)]
        for i in range(0,4):
            if cleavage2[0] == binding2[0]:
                continue
            try:
                if cleavage2[i]==binding2[i]:
                    match_counts[4+i] += 1
                total_counts[4+i] += 1
            except:
                pass
        for i in range(0, 4):
            if cleavage3[-1] == binding1[-1]:
                print(protein,row[1]['peptide'],start1,end1,start2,end2)
                continue
            try:
                if cleavage3[-1-i] == binding1[-1-i]:
                    match_counts[3-i] += 1
                total_counts[3-i] += 1
            except:
                pass
        pairs.append(cleavage3[-1]+binding1[-1])
        cleavage_bind1_bind2.append(cleavage3[-1]+binding1[-1]+binding2[0])
        try:
            if len(binding2[0:3])==3:
                bind_right_triple.append(binding2[0:3])
        except:
            pass
    df = pd.Series(pairs)
    counts = df.value_counts().reset_index()
    counts.columns=['pair','count']
    #counts = counts[(counts['pair'] == 'CS')|(counts['pair']=='SC')]
    #counts = counts[counts['count']>10]
    counts['expected'] = counts['pair'].apply(lambda x: aa_dict[x[0]]*aa_dict_spliced_no_bind1[x[1]]*len(pairs))
    counts['multiple'] = counts['count'].divide(counts['expected'])
    for i in possible_pairs:
        if i not in list(counts['pair']):
            counts = counts.append(pd.DataFrame([[i,0,0,0]],columns = ['pair', 'count', 'expected', 'multiple']),ignore_index=True)
    counts['expected'] = counts['pair'].apply(lambda x: aa_dict[x[0]] * aa_dict_spliced_no_bind1[x[1]] * len(pairs))
    print(counts.sort_values(by='multiple', ascending=False))
    print(match_counts,total_counts)
    print([x/y for x,y in zip(match_counts,total_counts)])
    counts.sort_values(by='multiple', ascending=False).to_csv('cleavage_to_bind1_comparison.csv')
    df_triple = pd.Series(cleavage_bind1_bind2)
    counts = df_triple.value_counts().reset_index()
    counts.columns = ['cleavage_bind1_bind2', 'count']
    counts['expected'] = counts['cleavage_bind1_bind2'].apply(lambda x: aa_dict[x[0]] * aa_dict_spliced_no_bind1[x[1]] * aa_dict_spliced[x[2]] * len(cleavage_bind1_bind2))
    counts['multiple'] = counts['count'].divide(counts['expected'])
    for i in possible_triples:
        if i not in list(counts['cleavage_bind1_bind2']):
            counts = counts.append(pd.DataFrame([[i,0,0,0]],columns = ['cleavage_bind1_bind2', 'count', 'expected', 'multiple']),ignore_index=True)
    counts['expected'] = counts['cleavage_bind1_bind2'].apply(lambda x: aa_dict[x[0]] * aa_dict_spliced_no_bind1[x[1]] * aa_dict_spliced[x[2]] * len(cleavage_bind1_bind2))
    print(counts.sort_values(by=['multiple','expected'],ascending=False))
    #counts = counts[(counts['cleavage_bind1_bind2'].str.startswith('RG'))|(counts['cleavage_bind1_bind2'].str.startswith('GR'))]
    #counts = counts[counts['count'] > 2]
    counts.sort_values(by=['multiple', 'expected'], ascending=False).to_csv('triple_comparison.csv')
    # Add categories
    counts['cleavage_bind1_bind2'] = counts['cleavage_bind1_bind2'].map(lambda x: aa_cats[x[0]]+aa_cats[x[1]]+aa_cats[x[2]])
    counts = counts[['cleavage_bind1_bind2','count','expected']].groupby('cleavage_bind1_bind2').sum()
    counts['ratio'] = counts.apply(lambda x:x[1]/(x[0]+0.00001),axis=1)
    counts.sort_values(by='ratio').to_csv('cats.csv')
    # Look at binding triples
    
    df_triple = pd.Series(bind_right_triple)
    counts = df_triple.value_counts().reset_index()
    counts.columns = ['cleavage_bind1_bind2', 'count']
    counts['expected'] = counts['cleavage_bind1_bind2'].apply(lambda x: aa_dict_spliced_no_bind1[x[0]] * aa_dict_spliced_no_bind1[x[1]] * aa_dict_spliced_no_bind1[x[2]] * len(bind_right_triple))
    counts['multiple'] = counts['count'].divide(counts['expected'])
    for i in possible_triples:
        if i not in list(counts['cleavage_bind1_bind2']):
            counts = counts.append(pd.DataFrame([[i,0,0,0]],columns = ['cleavage_bind1_bind2', 'count', 'expected', 'multiple']),ignore_index=True)
    counts['expected'] = counts['cleavage_bind1_bind2'].apply(lambda x: aa_dict[x[0]] * aa_dict_spliced_no_bind1[x[1]] * aa_dict_spliced[x[2]] * len(cleavage_bind1_bind2))
    print(counts.sort_values(by=['multiple','expected'],ascending=False))
    #counts = counts[(counts['cleavage_bind1_bind2'].str.startswith('RG'))|(counts['cleavage_bind1_bind2'].str.startswith('GR'))]
    #counts = counts[counts['count'] > 2]
    counts.sort_values(by=['multiple', 'expected'], ascending=False).to_csv('triple_comparison_binding.csv')
            
        
def plot_bind1(cell_line=''):
    df = pd.DataFrame.from_csv('all_spliced_peptides_with_cell_lines.csv')
    if cell_line:
        df = df[df['cell_line']==cell_line]
    df.drop('cell_line',axis=1,inplace=True)
    df.drop_duplicates(subset='peptide',keep=False,inplace=True)
    nonambiguous = set(pd.DataFrame.from_csv('all_spliced_peptides_nonambiguous.csv')['0'])
    for item in df.peptide:
        if item not in nonambiguous:
            df = df[df.peptide!=item]
     
        

#netchop_cleavage_predictions_unspliced()
#plot_9mer_sequence2(cell_line='C1R')
plot_all_9mers('GRLCL')
#generate_binding_pairs('left')
#plot_single_residue_between_peptides()
#peptide_to_cleavage_comparison()