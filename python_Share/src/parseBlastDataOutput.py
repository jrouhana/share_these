'''
Created on Aug 6, 2017

@author: jrouhana
'''
import re

###NOTE THAT FILE WAS PULLING DATA OFF A SERVER, WHICH IT DOES NOT HAVE ACCESS TO HERE###

# Class for blast objects
class BLAST(object):
    def __init__(self, blast_line):
        # Regex for blast terms
        blast_matches = re.search(r"""^(?P<transcript_id>.*?)\|.*?\t
                                gi\|.*?\|sp\|(?P<swissprot_id>.*?)\..*?\|.*?\t
                                (?P<pident>.*?)\t""", blast_line, re.M | re.S | re.X)

        # Define attributes 
        self.transcript_id = blast_matches.group("transcript_id")
        self.swissprot_id = blast_matches.group("swissprot_id")
        self.pident = float(blast_matches.group("pident"))

# Function that tests whether BLAST object identity > 95    
def blast_identity(blast_object):    
    return(blast_object.pident > 95)
    
# Create list of BLAST objects
blast_file = open("/scratch/RNASeq/blastp.outfmt6")
list_blast_objects = [BLAST(line)for line in blast_file.readlines()]
blast_file.close

# Load blast into dictionary if pident > 95
blast_dictionary = {blast_object.transcript_id:blast_object.swissprot_id for blast_object in list_blast_objects if blast_identity(blast_object)}

# Class for DE matrix objects
class MATRIX(object):
    def __init__(self, matrix_line):
        matrix_matches = re.search(r"""^(?P<transcript_id>.*?)\t
                                    (?P<ds>.*?)\t
                                    (?P<hs>.*?)\t
                                    (?P<log>.*?)\t
                                    (?P<plat>.*?)$""", matrix_line, re.M | re.S | re.X) 
        # Define attributes
        self.protein = matrix_matches.group("transcript_id")
        #If transcript in dictionary, replace with protein
        if self.protein in blast_dictionary:
            self.protein = blast_dictionary.get(self.protein)
        self.ds = float(matrix_matches.group("ds"))
        self.hs = float(matrix_matches.group("hs"))
        self.log = float(matrix_matches.group("log"))
        self.plat = float(matrix_matches.group("plat"))

# Create a list of DE objects
de_file = open("/scratch/RNASeq/diffExpr.P1e-3_C2.matrix")
list_de_objects = [MATRIX(line) for line in de_file.readlines()]
de_file.close()

# Output file
write_file = open("dataOutput.txt", "w")

# Print output
for de_object in list_de_objects:
        write_file.write(de_object.protein + "\t" + de_object.ds + "\t" + de_object.hs + "\t" + de_object.log + "\t" + de_object.plat + "\n")
write_file.close()
