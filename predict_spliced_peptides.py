from keras.models import load_model
from Bio import SeqIO
import pandas as pd

def predict_peptides(fasta):
    model = load_model('neural_net.h5')
    peptide_list = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        protein_seq = seq_record.seq
        protein_name = seq_record.id
        length = 9
        for i in range(len(protein_seq)):
            for j in range(1,9):
                length1 = j
                length2 = length-j
                sequence1 = protein_seq[i:i+j]
                cleavage1 = protein_seq[i-3,i+2]
                cleavage2 = protein_seq[i+j - 2, i+j + 3]
                for k in range(0,20):
                    distance = k
                    if i+9+k<len(protein_seq):
                        sequence2 = protein_seq[i+j+k:i+9+k]
                        cleavage3 = protein_seq[i+k - 3, i+k + 2]
                        cleavage4 = protein_seq[i + j +k - 2, i + j + 3]
                        for reversed in (True,False):
                            if reversed == True:
                                length1, length2, sequence1, sequence2, cleavage1, cleavage2, cleavage3, cleavage4 \
                                    = [length2, length1, sequence2, sequence1, cleavage3, cleavage4, cleavage1, cleavage2]
                            sequence = sequence1 + sequence2
                            peptide_list.append([sequence,length1,length2,l])
        df = pd.DataFrame()
                    
        
    
    
        model.predict_proba()
