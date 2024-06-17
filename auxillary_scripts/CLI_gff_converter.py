### 9/7/23
### Jessica L Albert
### CLI file - HIV Isoform Filtering
### last updated 9/12/23 by JLA

import argparse
import os
from gff_converter import *
#  
def main():
    # create parser object
    parser = argparse.ArgumentParser(prog = "GFF_converter",
                                     formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description =('''GFF Converter
Author: Jessica L Albert'''))
 
    # defining arguments for parser object
    parser.add_argument("input", type = str, 
                        metavar = "input_file_name", default = None,
                        help = "Designates input file to be filtered. This is required.")
    
    parser.add_argument("output", type = str, 
                        metavar = "output_file_prefix", default = None,
                        help = "Designates output file prefix. This is required.")
    
    parser.add_argument("ref", type = str,
                        metavar = "ref_file", default = None,
                        help = "Designates reference file name. This should be a fasta file. This is required.")
     
    
    # parse the arguments from standard input
    args = parser.parse_args()

    input_file = args.input
    
    
    output_file = args.output

    ref_file = args.ref
     
    
    convert_gff(input_file, output_file, ref_file)
 
if __name__ == "__main__":
    # calling the main function
    main()






