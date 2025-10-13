import csv
import sys
import os

################################ MAXIMIZATION LIST #############################

MAXIMIZATION_LIST = []

#############################################################################

#------------------------- BEGINNING OF FUNCTIONS --------------------------

### The function that will read the fimo.tsv file cleanly ###

def load_rows(path):
    with open(path, newline='') as fh:
        rdr = csv.reader(fh, delimiter='\t')   # ← TAB, not comma
        next(rdr, None)                        # skip header safely
        # keep only rows that have at least 3 cols (motif, alt-id, sequence)
        return [row for row in rdr if len(row) >= 3]

############################################################################

###  The function that will build the dictionary of motifs ###
def dictionaryBuilder(rows):
    motifDict = {}
    for row in rows:
        motifID = row[1]
        sequence = row[2]
        if motifID in motifDict:
            motifDict[motifID].append(sequence)
        else:
            motifDict[motifID] = [sequence]
    return motifDict

### An option with the dictionary built is to print it to another file,
### but what I advise is to simply gather the most common motif by 
### sequence for maximization analysis. 
############################################################################

### The function to write top motifs to a CSV file and their corresponding
### sequences in the following format:
### List_of_motifs = [motif1,motif2,motif3]
### >Motif 1, {number of sequences it appears in}
### [Sequence 1, Sequence 2, Sequence 3, Sequence 4]
### >Motif 2, {number of sequences it appears in}
### [Sequence 1, Sequence 2, Sequence 3]
### The goal is for this file to be intuitive ###

def fileBuilder(new_file_name,top_motif,top_value): 
    if os.path.exists(new_file_name):
        with open(new_file_name, "a") as file:
            length =  len(top_value)
            header = f">{top_motif}, number of sequences: {length}\n"
            file.write(header)
            file.write(str(top_value)+"\n")        
    else:
        with open(new_file_name, "w") as file:
            length =  len(top_value)
            header = f">{top_motif}, number of sequences: {length}\n"
            file.write(header)
            file.write(str(top_value)+"\n")


# -------------------------- END OF FUNCTIONS ------------------------------


# ------------------------- BEGINNING OF CALLS -----------------------------

### Loading the fimo.tsv file needed for the maximization analysis ###
if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} <fimo.tsv>")

rows = load_rows(sys.argv[1])
MAXIMIZATION_LIST_file = sys.argv[2]
COUNT_OF_RUNS = int(sys.argv[3])
############################################################################

### Building the motif and corresponding sequences file ###

first_draft_dictionary = dictionaryBuilder(rows)

################## Populating the maximization list ########################

with open(MAXIMIZATION_LIST_file, "r") as file:
    for line in file:
        motif = line.strip()
        if motif:  # Ensure the line is not empty
            MAXIMIZATION_LIST.append(motif)

############################################################################

### Extracting the top motif from the fimo.tsv output. This includes,
### the common motif name, the number of sequences the motif appears in,
### and a list of sequences the motif appears in For maximization anaylysis
### the practice will be to omit the motifs that have already been gathered,
### thus the top motif will not be the most common motif but rather, the 2nd,
### 3rd and so on from the top based on the number of loops through the
### pipeline


motif_to_add = None
value_to_add = []

for key, value in first_draft_dictionary.items():
    if len(value) > len(value_to_add) and key not in MAXIMIZATION_LIST:
        motif_to_add = key
        value_to_add = value

# Deduplicate while preserving order
seen = set()
unique_values = []
for seq in value_to_add:
    if seq not in seen:
        seen.add(seq)
        unique_values.append(seq)

value_to_add = unique_values


print("Most common motif:", motif_to_add)
print("Appears in", len(value_to_add), "sequences")
print("Sequences:", value_to_add)


### Building the file that will contain motifs and sequences, the one argument
### change is the name of the file by sample. Also building a file that will
### contain only the sequences of interest from the list of values


COUNT_OF_RUNS = str(COUNT_OF_RUNS)

with open(f"Test_values_{COUNT_OF_RUNS}.txt", "w", newline="") as file:
    writer = csv.writer(file)
    for row in value_to_add:
        writer.writerow([row])

fileBuilder(f"Test_counts_{COUNT_OF_RUNS}.txt",motif_to_add,value_to_add)

############################################################################




