import csv
import pandas as pd
'''
Add columns with netchop scores to spliced peptide data
'''

def new_columns(input,output):
    df = pd.DataFrame()
    with open(input,'r') as infile, open('netchopoutput.csv','r') as f, open(output,'w') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        for row in reader:
            protein = 'sp|' + row[1] + '|'
            [start1,end1,start2,end2] = row[2:6]
            reader2 = csv.reader(f)
            for row2 in reader2:
                if (row2[4] == protein) and (row2[0]==start1):
                    print(row)
                    
def new_columns_df(input,output):
    netchop_df = pd.read_csv('netchopoutput.csv')
    netchop_df.columns=['position','AA','S','probability','protein']
    output_df = pd.DataFrame()
    with open(input, 'r') as infile:
        reader = csv.reader(infile)
        for row in reader:
            protein = 'sp|' + row[1] + '|'
            data = row[2:6]
            [start1,end1,start2,end2] = [int(i) for i in data]
            protein_df = netchop_df[(netchop_df.protein == protein)]
            peptide1_df = protein_df[(protein_df.position >= start1) & (protein_df.position <= end1)]
            peptide2_df = protein_df[(protein_df.position >= start2) & (protein_df.position <= end2)]
            print(output_df)
        print(output_df)
            
new_columns_df('csv_files_parsed_peptides/GRLCL_spliced_peptides_new.csv', 'GRLCL_2D_peptides_with_cleavage_values.csv')
        