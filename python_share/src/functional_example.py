'''
Created on Jul 31, 2017
Updated on June 21, 2018

@author: jrouhana
'''
from __init__ import Prot_Dictionary, Tuple_Create, Tuple_To_Tab
import os

# Make sure input file exists
blast_file = "../data/blastp.outfmt6"
if not os.path.isfile(blast_file):
    print("The example blast data file is missing.")
    exit()

de_file = "../data/diffExpr.P1e-3_C2.matrix"
if not os.path.isfile(de_file):
    print("The example DE data file is missing.")
    exit()

# Create output directory if it does not exist
out_dir = "../output/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)


#Open Blast file
Blast_File = open(blast_file)
Swiss_Prot_Dictionary = Prot_Dictionary(Blast_File)
# close blast read file
Blast_File.close()

# Create a list of DE lines
DE_File = open(de_file)
List_DE_Lines = [Tuple_Create(line, Swiss_Prot_Dictionary) for line in DE_File.readlines()]
DE_File.close()




# Creates list of tabbed DE lines
List_Tabbed_DE = map(Tuple_To_Tab, List_DE_Lines)

# Print to output file, merging lines on \n, and tuples on \t
Output = open(out_dir + "swissProtStressData.txt", "w")
Output.write("\n".join(map(str, List_Tabbed_DE))+"\n")
Output.close()
