import pandas as pd
import matplotlib.pyplot as plt

def proline_analysis():
    df = pd.DataFrame.from_csv('GRLCL_spliced_9mers_binding_only_best.csv')
    df[df['peptide'].str[0]=='P'].boxplot('ic50',showfliers=False)
    plt.show()
    df[df['peptide'].str[0] != 'P'].boxplot('ic50',showfliers=False)
    plt.show()
    print(df[df['peptide'].str[0] != 'P']['ic50'].median())
    print(df[df['peptide'].str[0]=='P']['ic50'].quantile([0.,0.975]))
    print(df[df['peptide'].str[0] != 'P']['ic50'].quantile([0.025, 0.975]))

def MHC_counts():
    df1 = pd.DataFrame.from_csv('GRLCL_spliced_9mers_binding_only_best.csv')
    print(df1['allele'].value_counts())
    
MHC_counts()
    