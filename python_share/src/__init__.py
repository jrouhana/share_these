import re

#Some useful classes

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








#Some useful functions

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

# Function to return tuple from DE line
def Tuple_Create (DE_line, Swiss_Prot_Dictionary):
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



