'''
Created on Jul 17, 2017
Updated June 21, 2018

@author: jrouhana
'''
# Module for handling Regex
import os
import re
from __init__ import GOTermRecord


# Make sure input file exists
go_file = "../data/go-basic.obo"
if not os.path.isfile(go_file):
    print("The example go data file is missing.")
    exit()

de_file = "../data/diffExpr.P1e-3_C2.matrix"
if not os.path.isfile(de_file):
    print("The example DE data file is missing.")
    exit()

# Create output directory if it does not exist
out_dir = "../output/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)


# Splits GO file into records to pass into 
# class builder. Defaults to "go-basic.obo"
def Parse_GO_File(GO_File=go_file):
    # Defines dictionary for return storage
    GO_Dictionary = {}
    # reads in whole file, then splits on pattern
    Complete_GO_File = open(GO_File)
    Complete_GO_Input = Complete_GO_File.read()
    # return a list of matches
    GO_Records = re.findall(r"^\[Term\].*?\n\n", Complete_GO_Input, re.M | re.S)
    Complete_GO_File.close()
    # Pass to class constructor, then print return
    for GO_Record in GO_Records:
        Record = GOTermRecord(GO_Record)
        (GO_Key, GO_Value) = Record.Return_All()
        # Fill dictionary
        GO_Dictionary[GO_Key] = GO_Value
    return (GO_Dictionary)

# Function to print output
def Print_GO_File(output_file=out_dir+"parsedGO.txt"):
    # Open output file for writing
    GO_Output = open(output_file, "w")
    GO_Dictionary = Parse_GO_File()
    for GO_Term in sorted(GO_Dictionary):
        GO_Output.write(str(GO_Dictionary[GO_Term]))
    GO_Output.close()
   
Print_GO_File()
