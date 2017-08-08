import matplotlib.pyplot as plt
import pandas as pd
aa_dict = {'S': 0.083207136799280254, 'V': 0.05969329362753354, 'A': 0.070178284858173401, 'C': 0.022970256771492441, 'Q': 0.047658794128419044, 'P': 0.063103624877869691, 'N': 0.035881471493242882, 'T': 0.053571477460770595, 'G': 0.065702594916323465, 'M': 0.021322187589502981, 'X': 3.3560167684260995e-06, 'L': 0.099593771833428171, 'W': 0.012186050151078161, 'E': 0.071012696606281039, 'F': 0.03651443391901419, 'Y': 0.026648097884764452, 'H': 0.026268603041240058, 'I': 0.043406014352977189, 'K': 0.057270072888268402, 'D': 0.047394198701097866, 'R': 0.056413582082473761}

def all_positions_group():
    df1 = pd.DataFrame.from_csv('splicing_rules/ratio_of_spliced_to_non_spliced_in_each_position.csv',index_col=0,header=0)
    aa_groups = {'S': 'ST', 'V':'ILV' , 'A':'AG', 'C':'C' , 'Q': 'NQ', 'P':'P' , 'N':'NQ', 'T': 'ST', 'G':'AG' , 'M':'M', 'L': 'ILV', 'W': 'FYW', 'E': 'E', 'F':'FYW', 'Y':'FYW' , 'H':'H' , 'I': 'ILV', 'K':'K' , 'D':'D' , 'R':'R'}
    f, axes = plt.subplots(3, 3, sharey=True)
    for i in range(3):
        for j in range(3):
            k = 3 * i + j + 1
            df2 = df1[(df1['pos'] == k)].sort_values(by='aa', ascending=True)[
                ['aa', 'observed_x', 'observed_y']].reset_index(drop=True)
            df2.columns = ['AminoAcid', 'Splice Frequency', 'Nonspliced Frequency']
            df2['Type'] = [aa_groups[x] for x in df2['AminoAcid']]
            df2 = df2.groupby('Type').sum()
            df2.reset_index(inplace=True)
            #print(df2)
            axe = axes[i, j]
            ax = df2.plot(ax=axe, kind='bar')
            ax.set_xticks(df2.index)
            ax.set_xticklabels(df2.Type)
    plt.show()
    
def all_positions():
    df1 = pd.DataFrame.from_csv('splicing_rules/ratio_of_spliced_to_non_spliced_in_each_position.csv', index_col=0,
                                header=0)
    df_aa = pd.DataFrame.from_dict(aa_dict,orient='index')
    df_aa.index.name = 'aa'
    df_aa.reset_index(inplace=True)
    df_aa.columns = ['aa','Proteome Frequency']
    f, axes = plt.subplots(3, 3, sharey=True)
    for i in range(3):
        for j in range(3):
            k = 3 * i + j + 1
            df2 = df1[(df1['pos'] == k)].sort_values(by='aa', ascending=True)[
                ['aa', 'observed_x', 'observed_y']].reset_index(drop=True)
            df2 = pd.merge(df2,df_aa, on='aa')
            print(df2)
            df2.columns = ['AminoAcid', 'Splice Frequency', 'Nonspliced Frequency','Proteome Frequency']
            axe = axes[i, j]
            ax = df2.plot(ax=axe, kind='bar')
            ax.set_xticks(df2.index)
            ax.set_xticklabels(df2.AminoAcid)
    plt.show()

def cleavage_positions():
    positions = ['left_left', 'left', 'central', 'right', 'right_right']
    titles = ['P3','P2','P1','P1\'','P2\'']
    f, axes = plt.subplots(5, 1, sharey=True)
    for i in range(5):
        pos = positions[i]
        df1 = pd.DataFrame.from_csv('splicing_rules/cleavage_frequencies_'+pos+'_residue_first_peptide.csv', index_col=None,
                                    header=0)
        df2 = pd.DataFrame.from_csv('splicing_rules/cleavage_frequencies_'+pos+'_residue_nonspliced.csv', index_col=None,
                                    header=0)
        df_aa = pd.DataFrame.from_dict(aa_dict, orient='index')
        df_aa.index.name = 'aa'
        df_aa.reset_index(inplace=True)
        df_aa.columns = ['aa', 'Proteome Frequency']
        df = pd.merge(df1,df2, on='aa')
        print(df)
        df = df[['aa','observed_x','observed_y']]
        df = pd.merge(df,df_aa, on='aa')
        df.columns = ['AminoAcid', 'Splice Frequency', 'Nonspliced Frequency','Proteome Frequency']
        df=df.sort_values(by='AminoAcid')
        axe = axes[i]
        ax = df.plot(ax=axe, kind='bar')
        ax.set_title(titles[i]+' residue')
        ax.set_xticks(df_aa.index)
        ax.set_xticklabels(df.AminoAcid)
    f.subplots_adjust(hspace=.5)
    plt.show()
    
def proline_graph():
    df = pd.DataFrame.from_csv('training_data_3_positive.csv',index_col=0)
    df = df[df['seq_1']=='P']
    df1 = pd.DataFrame.from_csv('splicing_rules/ratio_of_spliced_to_non_spliced_in_each_position.csv', index_col=0,
                                header=0)
    df_aa = pd.DataFrame.from_dict(aa_dict, orient='index')
    df_aa.index.name = 'aa'
    df_aa.reset_index(inplace=True)
    df_aa.columns = ['aa', 'Proteome Frequency']
    f, axes = plt.subplots(8, sharey=True, sharex=True)
    for i in range(2,10):
        df2 = df1[(df1['pos'] == i)].sort_values(by='aa', ascending=True)[
            ['aa', 'observed_y']].reset_index(drop=True)
        df3 = df['seq_'+str(i)].value_counts()
        total = df3.sum()
        df3 = df3/total
        df3 = df3.to_frame()
        df3 = df3.reset_index()
        df3.columns = ['aa','prolinefreq']
        df2 = pd.merge(df2, df_aa, on='aa')
        df2 = pd.merge(df2,df3,on='aa')
        print(df2)
        df2.columns = ['AminoAcid', 'Nonspliced Frequency', 'Proteome Frequency','Proline Frequency']
        axe = axes[i-2]
        ax = df2.plot(ax=axe, kind='bar')
        ax.set_xticks(df2.index)
        ax.set_xticklabels(df2.AminoAcid)
    plt.show()

def plot_counts():
    df_bind1 = pd.DataFrame.from_csv('splicing_rules/bind_frequencies_left_residue_latest.csv')['observed'].to_frame()
    df_bind2 = pd.DataFrame.from_csv('splicing_rules/bind_frequencies_right_residue_latest.csv')['observed'].to_frame()
    df_nonspliced = pd.DataFrame.from_csv('splicing_rules/residue_count_nonspliced.csv',header=None)
    df_nonspliced.columns = ['count']
    df_nonspliced['frequency'] = df_nonspliced/df_nonspliced.sum()
    df = pd.merge(df_nonspliced['frequency'].to_frame(),df_bind1,left_index=True,right_index=True)
    df.sort_index(inplace=True)
    df.plot(kind='bar')
    plt.show()

cleavage_positions()