#!/bin/bash -l

#SBATCH --partition=general-compute 
#SBATCH --qos=general-compute 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --tasks-per-node=1 
#SBATCH --constraint="[ICE-LAKE-IB|CASCADE-LAKE-IB]"
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --mail-user=gentigas@buffalo.edu
#SBATCH --mail-type=END
#SBATCH --output=logs/fimo_%x_%j.out
#SBATCH --error=logs/fimo_%x_%j.err

### The modules needed to run FIMO

module load gcc/11.2.0 openmpi/4.1.1 ucx/1.11.2 meme

##############################################################
### INPUT FILES in order: test_sequence_fasta, the_database, the_number_of_runs
### The variables for analysis and naming

Maximization_Run=0
testSeq=$1
dataBase=$2
baseName=$(basename "$testSeq" .fa)
theDate=$(date +%Y-%m-%d)
outdir="${baseName}.${theDate}.${Maximization_Run}.output"
Maximization_List=()
Selected_motifs="./Selected_motifs.txt"
: > "$Selected_motifs"
Selected_motifs_report="./Selected_motifs_report.txt"
Count_Script="./Motif_Sequence_Counter_12_6_2025.py"
Number_of_Pipeline_Runs=$3
Number_of_sequences=$(grep -c "^>" "$testSeq")
### Optional
Common_Motifs_File=$4
echo "Running FIMO on $testSeq" | tee -a "${baseName}.run.log"

##############################################################

########## Finding number of sequences in fasta ##############
## The purpose of this number is to decide whether the merging
## version of this pipeline will be run, or the regular FIMO
## friendly version where Fasta is less than 1000 sequences

NumberOfSeqs=$(grep -c "^>" "$testSeq")

#----------------------------------------------------------------

### Numbering each sequence in the Fasta file ran through FIMO
# Only number if not already numbered
if grep -q "^>[^_]*_[0-9]\+:" "$testSeq"; then
    echo "Sequences already numbered, skipping numbering step"
else
    echo "Numbering sequences in $testSeq"
    tmp=$(mktemp)
    awk '
        BEGIN { OFS="" }          # keep formatting
        /^>/ {                    # header line?
            sub(/^>/,"",$0)       # drop the leading ">"
            split($0,a,":")       # a[1] = chromosome, a[2] = coords
            ++idx[a[1]]           # bump per-chromosome counter
            $0=">"a[1]"_"idx[a[1]]":"a[2]
        }
        { print }
    ' "$testSeq" > "$tmp" && mv "$tmp" "$testSeq"
fi

##############################################################

#----------------------------------------------------------------

### Beginning of Functions

### The FIMO batch submission function
runFimo() {
	fimo --o "${outdir}" \
    	"${dataBase}" \
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

### Splitting function for FASTA file if it is too large for single FIMO run

split_fasta_by_count() {
  local infile=$1
  local base=$(basename "$infile" .fa)
  awk -v base="$base" '
    BEGIN { seq=0; out="" }
    /^>/ {
      seq++
      if (seq % 400 == 1) {
        idx = int((seq - 1) / 400) + 1
        out = sprintf("%s_%d.fa", base, idx)
        close(out)
      }
    }
    { print >> out }
  ' "$infile"
}

### Merging function for fimo output used if the FASTA file is too large for FIMO

merge_fimo_tsv_files() {
  local outfile=$1
  shift
  awk '
    FNR == 1 {
      if (!header_printed) {
        print
        header_printed = 1
      }
      next
    }
    /^#/ { next }
    {
      if (!seen[$0]++) print
    }
  ' "$@" > "$outfile"
}

### End of Functions

#---------------------------------------------------------------------

### Beginning of function calls

### First condition: Is the number of sequences in the FASTA file lower
### than 800, an arbitrary number that I find works on FIMO most of the 
### time. If true, proceeds with maximization analysis

if [ "$NumberOfSeqs" -lt 400 ]; then
	echo "Processing Fasta file that has less than 400 sequences"
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
				"$Common_Motifs_File"\
				> "${baseName}_${Maximization_Run}_motif_counts"
			    Count_File="${baseName}_${Maximization_Run}_motif_counts"
			    echo "Counting complete. Output: $Count_File" | tee -a "${baseName}.run.log"
			    line=$(grep "Most common motif:" "$Count_File")
			    motif=$(echo "$line" | awk '{print $4}')
			    line2=$(grep "Appears in" "$Count_File")
			    motif_count=$(echo "$line2" | awk '{print $3}')
			    motif_percentage=$(awk "BEGIN {printf \"%.2f\", ($motif_count/$Number_of_sequences)*100}")
			    Maximization_List+=("$motif")
			    echo "$motif" >> "$Selected_motifs"
			    echo "$motif appears in ($motif_percentage%) of sequences" >> "$Selected_motifs_report" 
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
			    "$Common_Motifs_File"\
			    > "${baseName}_${Maximization_Run}_motif_counts"
			    Count_File="${baseName}_${Maximization_Run}_motif_counts"
			    echo "Counting complete. Output: $Count_File" | tee -a "${baseName}.run.log"
			    line=$(grep "Most common motif:" "$Count_File")
			    motif=$(echo "$line" | awk '{print $4}')
			    line2=$(grep "Appears in" "$Count_File")
			    motif_count=$(echo "$line2" | awk '{print $3}')
			    motif_percentage=$(awk "BEGIN {printf \"%.2f\", ($motif_count/$Number_of_sequences)*100}")
			    Maximization_List+=("$motif")
			    echo "$motif" >> "$Selected_motifs"
			    echo "$motif appears in ($motif_percentage%) of sequences" >> "$Selected_motifs_report"
			    echo "Current Maximization List = ${Maximization_List[@]}" | tee -a "${baseName}.run.${Maximization_Run}.log"
			    echo "This concluded Maximization Run ${Maximization_Run}" | tee -a "${baseName}.run.${Maximization_Run}.log"
			    Maximization_Run=$((Maximization_Run + 1))
			else
			    echo "ERROR: ${outdir}/fimo.tsv not found" >&2 | tee -a "${baseName}.run.${Maximization_Run}.log"
			    exit 1
			fi
		fi
	done

### Second condition: If the number of sequences in the FASTA file is greater
### than 800, then file is first split in two a number of files, 800 sequences
### in length each. This allows the same file to be run through FIMO in chunks. 

elif [ "$NumberOfSeqs" -gt 400 ]; then
	echo "Processing Fasta file that has more than 400 sequences"

	# actually split the fasta
	split_fasta_by_count "$testSeq"
	base=$(basename "$testSeq" .fa)
	shopt -s nullglob
	files=(${base}_*.fa)
	shopt -u nullglob

	if [ ${#files[@]} -eq 0 ]; then
		echo "ERROR: no chunk files found for ${base}_*.fa" >&2
		exit 1
	fi

	mkdir -p "$outdir"

	for split_fasta in "${files[@]}"; do
		chunk_tag=$(basename "$split_fasta" .fa)
		echo "Running FIMO on $split_fasta"
		chunk_outdir="${outdir}/${chunk_tag}"
		mkdir -p "$chunk_outdir"

		fimo --oc "$chunk_outdir" \
			"$dataBase" \
			"$split_fasta" \
			2> "logs/fimoTest.${chunk_tag}.${theDate}.log"
	done

	# Merge all chunk fimo.tsv files into the main output directory
	merge_fimo_tsv_files "${outdir}/fimo_merged.tsv" "${outdir}"/*/fimo.tsv
	
	# First count analysis (Maximization_Run = 0)
	if [ -s "${outdir}/fimo_merged.tsv" ]; then
		echo "Running Motif_Sequence_Counter.py on ${outdir}/fimo_merged.tsv"
		python "$Count_Script" \
			"${outdir}/fimo_merged.tsv" \
			"$Selected_motifs" \
			"$((Maximization_Run+1))" \
			"$Common_Motifs_File"\
			> "${baseName}_${Maximization_Run}_motif_counts"
		Count_File="${baseName}_${Maximization_Run}_motif_counts"
		echo "Counting complete. Output: $Count_File" | tee -a "${baseName}.run.log"
		line=$(grep "Most common motif:" "$Count_File")
		motif=$(echo "$line" | awk '{print $4}')
		line2=$(grep "Appears in" "$Count_File")
		motif_count=$(echo "$line2" | awk '{print $3}')
		motif_percentage=$(awk "BEGIN {printf \"%.2f\", ($motif_count/$Number_of_sequences)*100}")
		Maximization_List+=("$motif")
		echo "$motif" >> "$Selected_motifs"
		echo "$motif appears in ($motif_percentage%) of sequences" >> "$Selected_motifs_report"
		echo "Current Maximization List = ${Maximization_List[@]}" | tee -a "${baseName}.run.${Maximization_Run}.log"
		echo "This concluded Maximization Run ${Maximization_Run}" | tee -a "${baseName}.run.${Maximization_Run}.log"
		Maximization_Run=$((Maximization_Run + 1))
	else
		echo "ERROR: ${outdir}/fimo_merged.tsv not found or empty" >&2
		exit 1
	fi

	# Subsequent runs (Maximization_Run > 0)
	while [ "$Maximization_Run" -le "$Number_of_Pipeline_Runs" ]; do
		outdir="${baseName}.${theDate}.${Maximization_Run}.output"
		Sequences_To_Keep="./Test_values_${Maximization_Run}.txt"
		outputFasta=$(filterFasta "$Sequences_To_Keep" "$testSeq")
		testSeq="$outputFasta"
		
		NumOfSeqs=$(grep -c "^>" "$testSeq")
		
		if [ "$NumOfSeqs" -gt 400 ]; then
			echo "Run $Maximization_Run: Still >400 sequences ($NumOfSeqs), splitting..."
			
			# Clean up old chunk files to avoid conflicts
			base_filtered=$(basename "$testSeq" .fa)
			rm -f "${base_filtered}_*.fa"
			
			split_fasta_by_count "$testSeq"
			base=$(basename "$testSeq" .fa)
			shopt -s nullglob
			files=(${base}_*.fa)
			shopt -u nullglob

			if [ ${#files[@]} -eq 0 ]; then
				echo "ERROR: no chunk files found for ${base}_*.fa" >&2
				exit 1
			fi

			mkdir -p "$outdir"

			for split_fasta in "${files[@]}"; do
				chunk_tag=$(basename "$split_fasta" .fa)
				echo "Running FIMO on $split_fasta"
				chunk_outdir="${outdir}/${chunk_tag}"
				mkdir -p "$chunk_outdir"

				fimo --oc "$chunk_outdir" \
					"$dataBase" \
					"$split_fasta" \
					2> "logs/fimoTest.${chunk_tag}.${theDate}.log"
			done
			
			merge_fimo_tsv_files "${outdir}/fimo_merged.tsv" "${outdir}"/*/fimo.tsv

			if [ -s "${outdir}/fimo_merged.tsv" ]; then
				echo "Running Motif_Sequence_Counter.py on ${outdir}/fimo_merged.tsv"
				python "$Count_Script" \
					"${outdir}/fimo_merged.tsv" \
					"$Selected_motifs" \
					"$((Maximization_Run+1))" \
					"$Common_Motifs_File"\
					> "${baseName}_${Maximization_Run}_motif_counts"
				Count_File="${baseName}_${Maximization_Run}_motif_counts"
				echo "Counting complete. Output: $Count_File" | tee -a "${baseName}.run.log"
				line=$(grep "Most common motif:" "$Count_File")
				motif=$(echo "$line" | awk '{print $4}')
				line2=$(grep "Appears in" "$Count_File")
				motif_count=$(echo "$line2" | awk '{print $3}')
				motif_percentage=$(awk "BEGIN {printf \"%.2f\", ($motif_count/$Number_of_sequences)*100}")
				Maximization_List+=("$motif")
				echo "$motif" >> "$Selected_motifs"
				echo "$motif appears in ($motif_percentage%) of sequences" >> "$Selected_motifs_report"
				echo "Current Maximization List = ${Maximization_List[@]}" | tee -a "${baseName}.run.${Maximization_Run}.log"
				echo "This concluded Maximization Run ${Maximization_Run}" | tee -a "${baseName}.run.${Maximization_Run}.log"
				Maximization_Run=$((Maximization_Run + 1))
			else
				echo "ERROR: ${outdir}/fimo_merged.tsv not found" >&2 | tee -a "${baseName}.run.${Maximization_Run}.log"
				exit 1
			fi
			
		else
			echo "Run $Maximization_Run: <400 sequences ($NumOfSeqs), running normally..."
			
			runFimo
			
			if [ -f "${outdir}/fimo.tsv" ]; then
				echo "Running Motif_Sequence_Counter.py on ${outdir}/fimo.tsv"
				python "$Count_Script" \
					"${outdir}/fimo.tsv" \
					"$Selected_motifs" \
					"$((Maximization_Run+1))" \
					"$Common_Motifs_File"\
					> "${baseName}_${Maximization_Run}_motif_counts"
				Count_File="${baseName}_${Maximization_Run}_motif_counts"
				echo "Counting complete. Output: $Count_File" | tee -a "${baseName}.run.log"
				line=$(grep "Most common motif:" "$Count_File")
				motif=$(echo "$line" | awk '{print $4}')
				line2=$(grep "Appears in" "$Count_File")
				motif_count=$(echo "$line2" | awk '{print $3}')
				motif_percentage=$(awk "BEGIN {printf \"%.2f\", ($motif_count/$Number_of_sequences)*100}")
				Maximization_List+=("$motif")
				echo "$motif" >> "$Selected_motifs"
				echo "$motif appears in ($motif_percentage%) of sequences" >> "$Selected_motifs_report"
				echo "Current Maximization List = ${Maximization_List[@]}" | tee -a "${baseName}.run.${Maximization_Run}.log"
				echo "This concluded Maximization Run ${Maximization_Run}" | tee -a "${baseName}.run.${Maximization_Run}.log"
				Maximization_Run=$((Maximization_Run + 1))
			else
				echo "ERROR: ${outdir}/fimo.tsv not found" >&2 | tee -a "${baseName}.run.${Maximization_Run}.log"
				exit 1
			fi
		fi
	done
fi
	
	
#----------------------------------------------------------------------------

### Cleaning directory
	
#### Moving all the results into new Directory ####
if [ "$Maximization_Run" -gt "$Number_of_Pipeline_Runs" ]; then
	mkdir "$baseName.MAXIMIZATION_RESULT"
	mkdir "$baseName.RUNS"
	resultDir="$baseName.MAXIMIZATION_RESULT"
	runs="$baseName.RUNS"
	logs="./logs/"
	mv *.txt "$resultDir"
	mv *_counts "$resultDir"
	mv *.log "$logs"
	mv *run* "$runs"
fi
