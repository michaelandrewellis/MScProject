import subprocess
import tempfile
import csv
import pandas as pd
netchop = "./netchop-3.1/netchop"
import time
import itertools
from Bio import SeqIO


class peptide:
    def __init__(self,sequence,start_position,end_position,length,protein):
        self.seq = sequence
        self.start = start_position
        self.end = end_position
        self.length = length
        self.protein = protein

class spliced_peptide:
    def __init__(self, peptide1,peptide2):
        self.seq = peptide1.seq + peptide2.seq
        self.seq1 = peptide1.seq
        self.start1 = peptide1.start
        self.end1 = peptide1.end
        self.length1 = peptide1.length
        self.seq2 = peptide2.seq
        self.start2 = peptide2.start
        self.end2 = peptide2.end
        self.length2 = peptide2.length
        self.length = self.length1 + self.length2
        self.protein = peptide1.protein

    
def parse_proteome():
    with open('netchopoutput20s.csv','w') as output_file:
        pass # delete contents of file
    #df = pd.DataFrame()
    #df.columns = ['protein','start','end','seq']
    with open('reviewedproteome.fasta') as f:
        text = f.read()
        fastas = ['>sp' + fasta for fasta in text.split('>sp')[1:]] # ignore first element - empty string
        count = 0
        time1 = time.time()
        for fasta in fastas:
            count += 1
            with tempfile.NamedTemporaryFile(suffix=".fsa", mode="w") as tmpfile:
                tmpfile.write(fasta)
                tmpfile.flush()
                output = subprocess.check_output([netchop, '-v 1', tmpfile.name]).decode()
                lines = output.splitlines()[17:-5]
                lines = [line.split() for line in lines] # list of lists containing
                df = pd.DataFrame(lines)
                with open("netchopoutput20s.csv",'a') as output_file:
                    df.to_csv(output_file,index=False,header=False)
                if count % 100 == 0:
                    time2 = time.time()
                    time_elapsed = time2-time1
                    avg_time = time_elapsed/count
                    print(count, time_elapsed, 'Estimated time remaining: ', avg_time*(20000-count),' seconds')
                #with open("testoutput.csv", "wb") as csv_file:
                 #   writer = csv.writer(csv_file)
                  #  writer.writerows(lines)
                #for line in lines:
                    # if line[] = add sequences to dataframe
                        # cleaved after S
                        
    
'''
def get_peptides(cleavage_table):
    """
    Returns table of peptides from table of NetChop output
    """
    # Maybe include cleavage probabilities
    protein = cleavage_table[0][4]
    start = 1
    seq = ''
    peptide_table = []
    for aa in cleavage_table:
        seq += aa[1]
        if aa[2] == 'S':
            end = int(aa[0])
            length = end-start+1
            peptide_table.append([seq,start,end,length,protein])
            seq = ''
            start = int(aa[0]) + 1
    return(peptide_table)
'''

def get_peptides(cleavage_table):
    """
    Returns table of peptides from table of NetChop output
    """
    # Maybe include cleavage probabilities
    protein = cleavage_table[0][4]
    start = 1
    seq = ''
    peptide_table = []
    for aa in cleavage_table:
        seq += aa[1]
        if aa[2] == 'S':
            end = int(aa[0])
            length = int(end)-int(start)+1
            pt = peptide(seq,start,end,length,protein)
            peptide_table.append(pt)
            seq = ''
            start = int(aa[0]) + 1
    return(peptide_table)

def get_spliced_peptides(peptide_table):
    spliced_peptide_list = []
    #allowable_lengths = [9]
    for tuple in itertools.permutations(peptide_table,2):
        peptide1 = tuple[0]
        peptide2 = tuple[1]
        if peptide1.length + peptide2.length ==9:
            sp = spliced_peptide(peptide1,peptide2)
            spliced_peptide_list.append(sp)
    return(spliced_peptide_list)
    
    
 
            
def csv_to_spliced_peptides(input,output):
    df = pd.read_csv(input, sep=',', header=None)
    proteinlist=[]
    with open('proteinlist.csv','r') as f:
        reader = csv.reader(f)
        for row in reader:
            proteinlist.extend(row)
            print(proteinlist)
    with open(output, 'w') as w:
        pass
    with open(output, 'a') as w:
        count =0
        for protein in proteinlist:
            cleavage_df = df.loc[df[4] == "sp|" + protein + "|"]
            if not cleavage_df.empty:
                cleavage_table = cleavage_df.values.tolist() # POSSIBLY SPLIT THIS INTO TWO FUNCTIONS
                peptide_table = get_peptides(cleavage_table)
                spliced_peptide_list = get_spliced_peptides(peptide_table)
                for peptide in spliced_peptide_list:
                    count += 1
                    if count % 100000 == 0:
                        print(count)
                    wr = csv.writer(w)
                    wr.writerow([peptide.seq,peptide.protein,peptide.start1,peptide.end1,peptide.start2,peptide.end2]) # add in all other columns
            
        
def unique_proteins(fasta,output):
    proteins = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        proteins.append(seq_record.id.split("|")[1])
        print(seq_record.id)
        print(seq_record.id.split("|"))
    with open(output,'w') as f:
        wr = csv.writer(f)
        wr.writerow(proteins)
        
# output = subprocess.check_output([netchop, "./netchop-3.1/test/test.fsa"]).decode()

def test_parsing():
    time1=time.time()
    with open("./netchop-3.1/test/test.fsa") as f:
        text = f.read()
        time2 =time.time()
        with tempfile.NamedTemporaryFile(suffix=".fsa", mode="w") as tmpfile:
            tmpfile.write(text)
            tmpfile.flush()
            output = subprocess.check_output([netchop, tmpfile.name]).decode()
            lines = output.splitlines()[17:-5]
            lines = [line.split() for line in lines]
            df = pd.DataFrame(lines)
            df.to_csv('smalltest.csv',header=False,index=False)
    csv_to_spliced_peptides('smalltest.csv','smalltestout.csv')

            
#csv_to_spliced_peptides('testoutput.csv','testsplicedoutput3.csv')
#unique_proteins('testoutput.csv','proteinlist.csv')
#unique_proteins('reviewedproteome.fasta','proteinlist.csv')
#parse_proteome()
