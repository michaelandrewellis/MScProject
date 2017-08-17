import glob
import re
import pandas as pd
import pyfaidx as pf
import random

def get_ambiguous_peptides():
    ambiguous_spliced_peptides = set()
    ambiguous_unspliced_peptides = set()
    for file in glob.glob('csv_files_original/*_spliced*'):
        with open(file,'r') as f:
            f.readline() # skip header
            lines = f.read().splitlines()
            lines = [re.split(',|;',x) for x in lines]
            ambiguous_peptides = [x[0] for x in lines if len(x)>2]
            ambiguous_spliced_peptides.update(set(ambiguous_peptides))
    for file in glob.glob('csv_files_original/*_unspliced*'):
        with open(file,'r') as f:
            f.readline() # skip header
            lines = f.read().splitlines()
            lines = [re.split(',|;',x) for x in lines]
            ambiguous_peptides = [x[0] for x in lines if len(x)>2]
            ambiguous_unspliced_peptides.update(set(ambiguous_peptides))
    return [ambiguous_spliced_peptides, ambiguous_unspliced_peptides]

def remove_ambiguous_peptides():
    '''
    Remove ambiguous peptides and save to file
    '''
    ambiguous_spliced_peptides, ambiguous_unspliced_peptides = get_ambiguous_peptides()
    df_spliced = pd.DataFrame.from_csv('all_spliced_peptides_with_cell_lines.csv')
    df_spliced= df_spliced[~df_spliced['peptide'].isin(ambiguous_spliced_peptides)]
    df_spliced.to_csv('all_spliced_peptides_unambiguous_with_cell_lines.csv')
    df_unspliced = pd.DataFrame.from_csv('all_unspliced_peptides_with_cell_lines.csv')
    df_unspliced = df_unspliced[~df_unspliced['peptide'].isin(ambiguous_unspliced_peptides)]
    df_unspliced.to_csv('all_unspliced_peptides_unambiguous_with_cell_lines.csv')
    
            
def create_invented_peptides(n,outputfile):
    proteins = pf.Fasta('../reviewedproteome.fasta', key_function=lambda x: x.split('|')[1])
    protein_names = [protein.name for protein in proteins]
    peptide_list = []
    for i in range(int(n)):
        passed = False
        while passed == False:
            protein = random.choice(protein_names)
            protein_length = proteins[protein].__len__()
            start1 = random.randint(1,protein_length)
            length1 = random.randint(1,8)
            length2 = 9-length1
            end1 = start1+length1-1
            distance = random.randint(0,19)
            reversed = random.choice([True,False])
            if reversed:
                end2 = start1-distance-1
                start2 = end2 - length2+1
            else:
                start2 = end1 + distance + 1
                end2 = start2 + length2 -1
            if (min(start1,start2)>=1) and (max(end1,end2) <= protein_length):
                if not ((distance==0)and(reversed==False)):
                    passed=True
                else:
                    pass
            else:
                pass
        sequence = proteins[protein][start1-1:end1].seq + proteins[protein][start2-1:end2].seq
        peptide_list.append([sequence,protein,start1,end1,start2,end2])
    df = pd.DataFrame(peptide_list)
    df.columns = ['peptide','protein','start1','end1','start2','end2']
    df.to_csv(outputfile)
            
def create_training_data_without_cleavage(fin,fout,cell_line=None):
    df_input = pd.DataFrame.from_csv(fin)
    df_input.drop_duplicates(inplace=True)
    if cell_line:
        df_input = df_input[df_input['cell_line']==cell_line]
    output = []
    for idx, row in df_input.iterrows():
        length1 = row.end1 - row.start1 +1
        length2 = row.end2 - row.start2 +1
        if length1+length2 != 9:
            continue
        binding1_1 = row.peptide[length1-1]
        if length1>1:
            binding1_2 = row.peptide[length1-2]
        else:
            binding1_2 = None
        binding2_1 = row.peptide[length1]
        if length2 > 1:
            binding2_2 = row.peptide[length1 + 1]
        else:
            binding2_2 = None
        reversed = row.start1>row.end2
        distance = min(abs(row.start2-row.end1-1),abs(row.start1-row.end2-1))
        output.append([row.peptide,length1,binding1_2,binding1_1,binding2_1,binding2_2,reversed,distance])
    df_output = pd.DataFrame(output,columns=['peptide','length1','binding1_2','binding1_1','binding2_1','binding2_2',
                                             'reversed','distance'])
    df_output[['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9']] = pd.DataFrame([list(x) for x in df_output['peptide']])
    df_output.drop('peptide',inplace=True,axis=1)
    df_output.to_csv(fout)
            
def create_training_data_with_cleavage(fin,fout,cell_line=None):
    proteins = pf.Fasta('../reviewedproteome.fasta', key_function=lambda x: x.split('|')[1])
    df_input = pd.DataFrame.from_csv(fin)
    if cell_line:
        df_input = df_input[df_input['cell_line']==cell_line]
    if 'cell_line' in df_input.columns:
        df_input.drop('cell_line',inplace=True,axis=1)
    df_input.drop_duplicates(inplace=True)
    output = []
    for idx, row in df_input.iterrows():
        length1 = row.end1 - row.start1 +1
        length2 = row.end2 - row.start2 +1
        if length1+length2 != 9:
            continue
        binding1_1 = row.peptide[length1-1]
        if length1>1:
            binding1_2 = row.peptide[length1-2]
        else:
            binding1_2 = None
        binding2_1 = row.peptide[length1]
        if length2 > 1:
            binding2_2 = row.peptide[length1 + 1]
        else:
            binding2_2 = None
        reversed = row.start1>row.end2
        distance = min(abs(row.start2-row.end1-1),abs(row.start1-row.end2-1))
        try:
            cleavage1 = proteins[row.protein][max(0, row.start1 - 3):row.start1 - 1].seq
            cleavage2 = proteins[row.protein][row.end1:min(row.end1 + 2, proteins[row.protein].__len__())].seq
            cleavage3 = proteins[row.protein][max(0, row.start2 - 3):row.start2 - 1].seq
            cleavage4 = proteins[row.protein][row.end2:min(row.end2 + 2, proteins[row.protein].__len__())].seq
        except:
            continue
        output.append([row.peptide, length1, binding1_2, binding1_1, binding2_1, binding2_2, reversed, distance,cleavage1,
                       cleavage2,cleavage3,cleavage4])
    df_output = pd.DataFrame(output,columns=['peptide','length1','binding1_2','binding1_1','binding2_1','binding2_2',
                                             'reversed','distance','cleavage1','cleavage2','cleavage3','cleavage4'])
    df_output[['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9']] = pd.DataFrame([list(x) for x in df_output['peptide']])
    for cleavage in ['cleavage1','cleavage2','cleavage3','cleavage4']:
        df_output[[cleavage+'_1',cleavage+'_2']] = pd.DataFrame([list(x) for x in df_output[cleavage]])
    df_output.drop('peptide',inplace=True,axis=1)
    df_output.drop(['cleavage1','cleavage2','cleavage3','cleavage4'],inplace=True,axis=1)
    df_output.to_csv(fout)

def add_ic50(fin,fout):
    df = pd.DataFrame.from_csv(fin)    
    df['peptide'] = df[['seq_1','seq_2','seq_3','seq_4','seq_5','seq_6','seq_7','seq_8','seq_9']].sum(axis=1)
    df_ic50 = pd.DataFrame.from_csv('GRLCL_spliced_9mers_binding_only_best.csv')
    print(df)

def create_allele_specific_training_data(allele,all_training_data_input_file,allele_training_data_output_file):
    df = pd.DataFrame.from_csv(all_training_data_input_file)
    df_binding = pd.DataFrame.from_csv('GRLCL_spliced_9mers_binding_only_best.csv')
    df_binding = df_binding[df_binding.ic50 < 500]
    df_binding = df_binding[df_binding.allele == allele]
    df['peptide'] = df[['seq_1', 'seq_2', 'seq_3', 'seq_4', 'seq_5', 'seq_6', 'seq_7', 'seq_8', 'seq_9']].sum(axis=1)
    df = df[df.peptide.isin(set(df_binding.peptide))]
    df.drop(['peptide'],inplace=True,axis=1)
    df.to_csv(allele_training_data_output_file)

def create_test_data_with_binding_above_500(all_training_data_input_file,output_file):
    df = pd.DataFrame.from_csv(all_training_data_input_file)
    df_binding = pd.DataFrame.from_csv('GRLCL_spliced_9mers_binding_only_best.csv')
    df_binding = df_binding[df_binding.ic50 > 500]
    df['peptide'] = df[['seq_1', 'seq_2', 'seq_3', 'seq_4', 'seq_5', 'seq_6', 'seq_7', 'seq_8', 'seq_9']].sum(axis=1)
    df = df[df.peptide.isin(set(df_binding.peptide))]
    df.drop(['peptide'], inplace=True, axis=1)
    df.to_csv(output_file)
    
def create_negative_training_data_with_binding(allele,all_training_data,output_file):
    df = pd.DataFrame.from_csv(all_training_data)
    df_binding = pd.DataFrame.from_csv('binding_predictions/invented_binding_'+allele+'_2.txt',sep='\t')
    df_binding = df_binding[df_binding.ic50 < 500]
    df['peptide'] = df[['seq_1', 'seq_2', 'seq_3', 'seq_4', 'seq_5', 'seq_6', 'seq_7', 'seq_8', 'seq_9']].sum(axis=1)
    df = df[df.peptide.isin(set(df_binding.peptide))]
    df.drop(['peptide'], inplace=True, axis=1)
    df.to_csv(output_file)

if __name__ == "__main__":
    #create_negative_training_data_with_binding('C0702', 'training_data/test_data_with_cleavage_negative.csv', 'training_data/test_data_negative_binders_C0702.csv')
    #create_invented_peptides(5000,'all_invented_peptides_2.csv')
    create_training_data_with_cleavage('all_spliced_peptides_unambiguous_with_cell_lines.csv','training_data/training_data_with_cleavage_positive_fibroblasts.csv','Fibroblasts')
    #add_ic50('training_data/training_data_with_cleavage_positive.csv','training_data/training_data_with_cleavage_and_ic50_positive.csv')
    #create_allele_specific_training_data('HLA-A*01:01','training_data/training_data_with_cleavage_positive_GRLCL.csv','training_data/training_data_with_cleavage_positive_A0101.csv')
    #create_test_data_with_binding_above_500('training_data/training_data_with_cleavage_positive_GRLCL.csv','training_data/test_data_positive_poor_binders.csv')
    #create_training_data_with_cleavage('all_invented_peptides_2.csv','training_data/test_data_with_cleavage_negative.csv')