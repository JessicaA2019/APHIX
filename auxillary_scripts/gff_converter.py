### 6/13/24
### Jessica L Albert
### GFF converter
### last updated 06/13/24 by JLA

import regex as re
import csv
import math
import pandas as pd
from Reference_aligner import *

def readFastqName(filename):
    #Reads FASTQ file and remove the special characters!
    
    with open(filename) as fh:
        name = fh.readline().rstrip() # read name
    return name

def convert_gff(input_file, output_file, ref_file):

    #caluculate reference genome info
    CIGAR = CIGAR_for_new_ref(ref_file)

    gff = pd.read_csv('auxillary_scripts/auxillary.gff', names=["ref", "User", "type", "start", "stop", "a", "b", "c", "name"], sep='\t')

    ref_name = readFastqName(ref_file)

    row = 0
    for ref in gff["ref"]:
        name = str(gff.iloc[row,0])
        if name == 'NL43_Complete':
            name = str(ref_name[1:])
            gff.iloc[row,0] = name
        elif name == '##NL43_Complete':
            name = "##"+str(ref_name[1:])
            gff.iloc[row,0] = name
        row = row + 1

    row = 0
    for number in gff.iloc[:,3]:    
        gff.iloc[row,3] = calc_shift(CIGAR, number)
        row = row + 1
        
    row = 0
    for number in gff.iloc[:,4]:
        gff.iloc[row,4] = calc_shift(CIGAR, number)
        row = row + 1

    #gff['start'].astype('Int64').astype('str')
    #gff['stop'].astype('Int64').astype('str')
    
    row = 0
    for entry in gff["start"]:
        name = str(gff.iloc[row,3])
        name = re.sub("\.(?:.)+$", "", name)
        if name != 'nan':
            gff.iloc[row,3] = name
        row = row + 1

    row = 0
    for entry in gff["stop"]:
        name = str(gff.iloc[row,4])
        name = re.sub("\.(?:.)+$", "", name)
        if name != 'nan':
            gff.iloc[row,4] = name
        row = row + 1


    gff.to_csv(output_file, index=False, sep='\t', header=False)


