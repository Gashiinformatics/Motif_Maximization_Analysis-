import csv
import sys
import os

################################ MAXIMIZATION LIST ############################

MAXIMIZATION_LIST = []

###############################################################################

################################### FLAG FOR LIST #############################
# This section populates the maximization list with previously seen common 
# motifs to filter them out.

# 1. Define the flag explicitly before using it
flag = True

if flag:
    # 2. Check if the 4th argument (the file path) was provided by Bash
    if len(sys.argv) > 4:
        optional_motif_filtering_file = sys.argv[4]

        # 3. Check if the file actually exists
        if os.path.exists(optional_motif_filtering_file):
            with open(optional_motif_filtering_file, "r") as f:
                motifs = f.readlines()
                for motif in motifs:
                    clean_motif = motif.strip()
                    # Only add if not empty
                    if clean_motif:
                        MAXIMIZATION_LIST.append(clean_motif)
        else:
            # It's okay if file doesn't exist, just print a warning
            print(f"Warning: File {optional_motif_filtering_file} does not exist.")
    else:
        # It's okay if the argument isn't there, just warn and continue
        print("Notice: Flag is True, but no common motifs file argument provided.")

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

############################################################################

### The function to write top motifs to a CSV file ###

def fileBuilder(new_file_name, top_motif, top_value): 
    # Use 'a' (append) if file exists, 'w' (write) if it doesn't
    mode = "a" if os.path.exists(new_file_name) else "w"
    
    with open(new_file_name, mode) as file:
        length = len(top_value)
        header = f">{top_motif}, number of sequences: {length}\n"
        file.write(header)
        file.write(str(top_value) + "\n")


# -------------------------- END OF FUNCTIONS ------------------------------


# ------------------------- BEGINNING OF CALLS -----------------------------

### Loading the fimo.tsv file needed for the maximization analysis ###

# Ensure we have the minimum required arguments
if len(sys.argv) < 4:
    sys.exit(f"Usage: {sys.argv[0]} <fimo.tsv> <exclusion_list> <run_number> [optional_common_motifs]")

rows = load_rows(sys.argv[1])
MAXIMIZATION_LIST_file = sys.argv[2]
COUNT_OF_RUNS = int(sys.argv[3])

############################################################################

### Building the motif and corresponding sequences file ###

first_draft_dictionary = dictionaryBuilder(rows)

################## Populating the maximization list from ARG 2 ##############
# This reads the accumulated list of motifs found IN THIS PIPELINE run
# (Separate from the optional file loaded at the top)

if os.path.exists(MAXIMIZATION_LIST_file):
    with open(MAXIMIZATION_LIST_file, "r") as file:
        for line in file:
            motif = line.strip()
            if motif:  # Ensure the line is not empty
                MAXIMIZATION_LIST.append(motif)

############################################################################

### Extracting the top motif from the fimo.tsv output. ###

motif_to_add = None
value_to_add = []

for key, value in first_draft_dictionary.items():
    # Only select if it has more sequences AND is not in our exclusion list
    if len(value) > len(value_to_add) and key not in MAXIMIZATION_LIST:
        motif_to_add = key
        value_to_add = value

# If we found nothing (e.g. everything was filtered out), handle gracefully
if motif_to_add is None:
    # We must print something so the Bash script's 'awk' doesn't crash
    print("Most common motif: None")
    print("Appears in 0 sequences")
    print("Sequences: []")
else:
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

    ### Building the output files ###

    COUNT_OF_RUNS = str(COUNT_OF_RUNS)

    with open(f"Test_values_{COUNT_OF_RUNS}.txt", "w", newline="") as file:
        writer = csv.writer(file)
        for row in value_to_add:
            writer.writerow([row])

    fileBuilder(f"Test_counts_{COUNT_OF_RUNS}.txt", motif_to_add, value_to_add)

############################################################################
