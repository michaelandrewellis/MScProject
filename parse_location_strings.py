import csv

def parse_location(input):
    with open(input,'r') as infile, open(input[:-4]+'_new.csv','w') as outfile:
        reader = csv.reader(infile)
        next(reader, None) # skip first line
        writer = csv.writer(outfile)
        for row in reader:
            row[1] = row[1].split(';')[0]
            protein = row[1].split('|')[1]
            locations = row[1].split(' ')[1].split('_')
            writer.writerow([row[0],protein]+locations)
        
parse_location('GRLCL_unspliced_peptides_1D_HCD.csv')