'''
Created on Jul 17, 2017
Updated June 21, 2018

@author: jrouhana
'''
# Module for handling Regex
import os
import re

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

# Define class for GO objects
class GOTermRecord(object):
    # Constructor
    def __init__(self, Record):
        # Regex for most terms
        GO_Matches = re.search(r"""^id:\s(?P<id>.*?)\n
                                ^name:\s(?P<name>.*?\n)
                                ^namespace:\s(?P<namespace>.*?\n)""", Record, re.M | re.S | re.X)
        
        # Regex for Isas matches
        GO_Isas_Matches = re.findall(r"^is_a:\s(?P<isa>.*?\n)", Record, re.M | re.S | re.X)
        
        # Define attributes 
        self.GO_Id = GO_Matches.group("id")
        GO_Name = GO_Matches.group("name")
        GO_Namespace = GO_Matches.group("namespace")
        if GO_Isas_Matches:
            # If the list is not empty, then
            GO_Isas_String = "\t".join(GO_Isas_Matches)
            self.Value_Seq = (GO_Namespace, GO_Name, GO_Isas_String)
        else: 
            self.Value_Seq = (GO_Namespace, GO_Name)

    # Method to return proper print output
    def Return_All(self):
        # Format for print
        Value = "\t".join(self.Value_Seq) 
        GO_Return = str(self.GO_Id + "\t" + Value + "\t\n")
        return(self.GO_Id, GO_Return)

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
