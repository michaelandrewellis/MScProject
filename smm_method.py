aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
import numpy as np
import os
from math import exp

allele = 'A-01:01'
length = 9

def allele_matrix(allele,length):
    with open(os.path.join(os.path.dirname(__file__), 'smm_matrix/HLA-' + allele +'-' + str(length)+'.txt'),'r') as f:
        f.readline() # skip header
        lines = f.readlines()
        matrix = []
        for i in range(0,20):
            matrix.append([float(j) for j in lines[i].split()[1:length+1]])
    array=np.array(matrix)
    return(array)

def peptide_matrix(peptide):
    matrix = []
    for letter in peptide:
        matrix.append(list(map(lambda x : int(x is letter),aa_list)))
    return(np.array(matrix))
        
    
def get_ic50(peptide,allele,length):
    allele_array = allele_matrix(allele,length)
    peptide_array = peptide_matrix(peptide)
    total = 0
    for i in range(0,length):
        total += np.dot(peptide_array[i,:],allele_array[:,i])
        print(np.dot(peptide_array[i,:],allele_array[:,i]))
    print(exp(total),total)
    
get_ic50('LSQQTLEEY','A-01:01',9)
