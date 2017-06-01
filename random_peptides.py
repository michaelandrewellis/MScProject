import pandas as pd
import random
'''
Generate 1500 random 9mers, 1000 10mers, 500 11mers, 250 12mers
'''
#length_dict = {9:1500,10:1000,11:500,12:250}
df = pd.read_csv('netchopoutput.csv')
df_length = df.shape[0]
table = []
length=9
number = 6000
#for length in length_dict:
for i in range(0,6000): #length_dict[length]
    length1 = random.randint(1, length-1)
    length2 = length - length1
    distance = random.randint(0, 25)
    reversed = random.choice([True, False])
    seq = ''
    cleavage1 = ''
    cleavage2 = ''
    cleavage3 = ''
    cleavage4 = ''
    r = random.randint(40, df_length-40)
    for j in range(r,r+length1):
        row = df.iloc[j]
        aa = row[1]
        seq += aa
    for j in range(r-3, r+2):
        row = df.iloc[j]
        aa = row[1]
        cleavage1 += aa
    for j in range(r+length1-3, r+length1+2):
        row = df.iloc[j]
        aa = row[1]
        cleavage2 += aa
    if reversed == True:
        r = r-distance - length2
    else:
        r = r + length1 + distance
    for j in range(r,r+length2):
        row = df.iloc[j]
        aa = row[1]
        seq += aa
    for j in range(r - 3, r + 2):
        row = df.iloc[j]
        aa = row[1]
        cleavage3 += aa
    for j in range(r + length2 - 3, r + length2 +2):
        row = df.iloc[j]
        aa = row[1]
        cleavage4 += aa
    #for j in range(0,9):
     #   r = random.randint(0,df_length-1)
      #  row = df.iloc[r]
       # aa = row[1]
        #seq += aa
    table.append([seq,cleavage1,cleavage2,cleavage3,cleavage4,distance,reversed,length1,length2,length1+length2])
df = pd.DataFrame(table)
df.columns = ['sequence','cleavage1','cleavage2','cleavage3','cleavage4','distance','reversed','length1','length2','length']
df.to_csv('6000_random_spliced_9mers_with_cleavage_sites.csv')
        