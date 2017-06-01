import pandas as pd
import numpy as np
aa_dict = {'S': 0.083207136799280254, 'V': 0.05969329362753354, 'A': 0.070178284858173401, 'C': 0.022970256771492441, 'Q': 0.047658794128419044, 'P': 0.063103624877869691, 'N': 0.035881471493242882, 'T': 0.053571477460770595, 'G': 0.065702594916323465, 'M': 0.021322187589502981, 'X': 3.3560167684260995e-06, 'L': 0.099593771833428171, 'W': 0.012186050151078161, 'E': 0.071012696606281039, 'F': 0.03651443391901419, 'Y': 0.026648097884764452, 'H': 0.026268603041240058, 'I': 0.043406014352977189, 'K': 0.057270072888268402, 'D': 0.047394198701097866, 'R': 0.056413582082473761}
amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

print(list(aa_dict.keys()))
print(list(aa_dict.values()))

seq_list=[]
for i in range(0,6000):
    seq = ''
    for i in range(0,9):
        aa = np.random.choice(list(aa_dict.keys()),1,p=list(aa_dict.values()))
        seq+=aa[0]
       
    seq_list.append(seq)