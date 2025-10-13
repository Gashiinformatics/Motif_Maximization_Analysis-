#!/bin/bash -l

#SBATCH --partition=general-compute 
#SBATCH --qos=general-compute 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --tasks-per-node=16 
#SBATCH --constraint="[ICE-LAKE-IB|CASCADE-LAKE-IB]"
#SBATCH --time=2:00:00
#SBATCH --mail-user=gentigas@buffalo.edu
#SBATCH --mail-type=END
#SBATCH --output=logs/fimo_%x_%j.out
#SBATCH --error=logs/fimo_%x_%j.err

### The modules needed to run FIMO

module load gcc/11.2.0 openmpi/4.1.1 ucx/1.11.2 meme

##############################################################

### The variables for analysis and naming

Maximization_Run=0
testSeq=$1
baseName=$(basename "$testSeq" .fa)
theDate=$(date +%Y-%m-%d)
outdir="${baseName}.${theDate}.${Maximization_Run}.output"
Maximization_List=()
Selected_motifs="./Selected_motifs.txt"
: > "$Selected_motifs"
Count_Script="./Motif_Sequence_Counter_VSCODE.py"
Number_of_Pipeline_Runs=4

echo "Running FIMO on $testSeq" | tee -a "${baseName}.run.log"

##############################################################

### The FIMO batch submission function

runFimo() {
	srun fimo --o "${outdir}" \
    	Drosophila_melanogaster.meme \
    	"${testSeq}" \
    	2> "fimoTest.${baseName}.${theDate}.log"
}

##############################################################

### Fasta Filter Function

filterFasta() {
    local keep_ids=$1
    local inputFasta=$2
    local outputFasta="${inputFasta%.*}_run_${Maximization_Run}.fa"
    awk -v keepfile="$keep_ids" '
        BEGIN {
            # Build unique keep-set; sanitize CRLF and blanks
            while ((getline k < keepfile) > 0) {
                sub(/\r$/, "", k);              # drop trailing CR (^M)
                if (k != "") keep[k]=1;         # dedup via assoc array
            }
        }
        /^>/ {
            id = $1;                            # first token of header
            sub(/^>/, "", id);                  # remove leading ">"
            sub(/:.*/, "", id);                 # drop suffix after first colon
            printseq = (id in keep);
        }
        printseq
    ' "$inputFasta" > "$outputFasta"
    echo "$outputFasta"
}


### The First Fimo submission

runFimo

###  The follow up motif counting analysis
while [ "$Maximization_Run" -le "$Number_of_Pipeline_Runs" ]
do
	if [ "$Maximization_Run" -eq 0 ]; then 
		# ---- Counting step right after FIMO ----

		if [ -f "${outdir}/fimo.tsv" ]; then
		    echo "Running Motif_Sequence_Counter.py on ${outdir}/fimo.tsv" 
		    python "$Count_Script" \
		    	"${outdir}/fimo.tsv" \
			"$Selected_motifs" \
			"$((Maximization_Run+1))"\
			 > "${baseName}_${Maximization_Run}_motif_counts"
	            Count_File="${baseName}_${Maximization_Run}_motif_counts"
		    line=$(grep "Most common motif:" "$Count_File")
                    motif=$(echo "$line" | awk '{print $4}')
                    Maximization_List+=("$motif")
                    echo "$motif" >> "$Selected_motifs"
		    echo "Current Maximization List = ${Maximization_List[@]}" | tee -a "${baseName}.run.${Maximization_Run}.log"
		    echo "This concluded Maximization Run ${Maximization_Run}" | tee -a "${baseName}.run.${Maximization_Run}.log"
			Maximization_Run=$((Maximization_Run + 1))
		else
		    echo "ERROR: ${outdir}/fimo.tsv not found" >&2 | tee -a "${baseName}.run.${Maximization_Run}.log"
		    exit 1
		fi
	elif [ "$Maximization_Run" -gt 0 ]; then
        	outdir="${baseName}.${theDate}.${Maximization_Run}.output"
		Sequences_To_Keep="./Test_values_${Maximization_Run}.txt"
    		outputFasta=$(filterFasta "$Sequences_To_Keep" "$testSeq")
    		testSeq="$outputFasta"
    		runFimo
		if [ -f "${outdir}/fimo.tsv" ]; then
		    echo "Running Motif_Sequence_Counter.py on ${outdir}/fimo.tsv"
		    python "$Count_Script" \
		    "${outdir}/fimo.tsv" \
		    "$Selected_motifs" \
		    "$((Maximization_Run+1))"\
		     > "${baseName}_${Maximization_Run}_motif_counts"
		    Count_File="${baseName}_${Maximization_Run}_motif_counts"
		    echo "Counting complete. Output: $Count_File" | tee -a "${baseName}.run.log"
		    line=$(grep "Most common motif:" "$Count_File")
		    motif=$(echo "$line" | awk '{print $4}')
		    Maximization_List+=("$motif")
		    echo "$motif" >> "$Selected_motifs"
		    echo "Current Maximization List = ${Maximization_List[@]}" | tee -a "${baseName}.run.${Maximization_Run}.log"
		    echo "This concluded Maximization Run ${Maximization_Run}" | tee -a "${baseName}.run.${Maximization_Run}.log"
		    Maximization_Run=$((Maximization_Run + 1))
                else
                    echo "ERROR: ${outdir}/fimo.tsv not found" >&2 | tee -a "${baseName}.run.${Maximization_Run}.log"
                    exit 1
		fi
	fi
done	
	
#### Moving all the results into new Directory ####
if [ "$Maximization_Run" -gt "$Number_of_Pipeline_Runs" ]; then
	mkdir "$baseName.MAXIMIZATION_RESULT"
	resultDir="$baseName.MAXIMIZATION_RESULT"
	logs="./logs/"
	mv *.txt "$resultDir"
	mv *_counts "$resultDir"
	mv *.log "$logs"
fi

