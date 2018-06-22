'''
Created on Jul 31, 2017

@author: jrouhana
'''

###NOTE THAT FILE WAS PULLING DATA FROM A SERVER, WHICH IT DOES NOT HAVE ACCESS TO HERE###

# Function to read Blast file and parse
def Parse_Blast (Blast_Line): 
    # split each line into a list of strings separated by tab
    Line_List = Blast_Line.rstrip("\n").split("\t")
    Column_One, Column_Two = Line_List[:2]
    # split column one by | to separate transcriptID
    Column_One_List = Column_One.split("|")
    # retrieve transcriptID from column_one
    TranscriptID = Column_One_List[0]
    # split column two by | to isolate what we need    
    Column_Two_List = Column_Two.split("|")
    SwissprotID_Version = Column_Two_List[3]
    # Remove version number
    SwissprotID_List = SwissprotID_Version.split(".")
    SwissprotID = SwissprotID_List[0]
    # Return values
    return(TranscriptID, SwissprotID)

# Function to create dictionary
def Prot_Dictionary (Blast_File):
    File = Blast_File
    # Define dictionary to store items in
    Swiss_Prot_Dictionary = {} 
    # Process BLAST file by line
    for Line in File:
        TranscriptID, Swiss_ProtID = Parse_Blast(Line)
        Swiss_Prot_Dictionary[TranscriptID] = Swiss_ProtID
    return Swiss_Prot_Dictionary


Blast_File = open("../data/blastp.outfmt6")
Swiss_Prot_Dictionary = Prot_Dictionary(Blast_File)
# close blast read file
Blast_File.close()

# Function to return tuple from DE line
def Tuple_Create (DE_line):
    # Splits on tabs, then returns tuple
    Line_List = DE_line.rstrip("\n").split("\t")
    if Line_List[0] in Swiss_Prot_Dictionary:
        Line_List[0] = Swiss_Prot_Dictionary.get(Line_List[0])
        return(tuple(Line_List))
    else:
        return(tuple(Line_List))

# Converts tuple to tab separated string
def Tuple_To_Tab (Tuple_Line):
    return "\t".join(Tuple_Line)
        
# Create a list of DE lines
DE_File = open("../data/diffExpr.P1e-3_C2.matrix")
List_DE_Lines = map(Tuple_Create, DE_File.readlines())
DE_File.close()

# Creates list of tabbed DE lines
List_Tabbed_DE = map(Tuple_To_Tab, List_DE_Lines)

# Print to output file, merging lines on \n, and tuples on \t
Output = open("../intermediate_files/swissProtStressData.txt", "w")
Output.write("\n".join(map(str, List_Tabbed_DE))+"\n")
Output.close()
