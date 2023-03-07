#!/bin/bash

## Data organisation
## create data structure in directories
## MiSeq data is different in raw seq data as genomescan data
## This wrapper is ONLY for MiSeq data.
## MiSeq data name structure "[SAMPLENAME]_S[NUMBER]_R[1,2]_[001 or empty].fastq.gz"
## This structure can be simplified to "[SAMPLENAME]_R[1,2].fastq.gz" for easy use

### manual input of some info regarding working directories
echo '2 options for analysis';
echo 'without cleanup of contaminant sequences ==> choose 1';
echo '                    removal of human, cat & dog sequences';
echo '                    with the use of a masked db to reduce false positives';
echo 'or';
echo 'with cleanup of contaminant sequences ==> choose 2';
echo '                    removal of human, cat & dog sequences';
echo '                    with the use of a masked db to reduce false positives';

SETTING="USER INPUT";
read -p "
Enter 1 for analysis without delete contaminant sequences, 
Enter 2 for analysis with delete contaminant sequences: " SETTING && [[ "$SETTING" == [1,2] ]] || exit 1


SETTING2="USER INPUT2";
read -p "
Enter --> 0 for k2_standard kraken2 database
Enter --> 1 for virus specific database
Enter --> 2 for k2_pluspfp kraken2 database: " SETTING2 && [[ "$SETTING2" == [0,1,2] ]] || exit 1


echo "$SETTING2";
selectDB="$SETTING2";

echo "$selectDB";



echo -e "Setting=$SETTING";

if [ "$SETTING" = 1 ];then

# settings used for the analysis
	echo 'Without cleanup of contaminants is used for analysis';
	echo 'alignment must have:'
	CLEANUP=0;

else
CONTAMINANT="USER INPUT";
read -p "
enter H for human 
enter B for bacteria 
enter C for cat
enter D for Dog
 " CONTAMINANT && [[ "$CONTAMINANT" == [H,h,B,b,C,c,D,d] ]] || exit 1
	CLEANUP=1;

# settings used for analysis
	echo 'With cleanup of contaminants is used for analysis';

fi

read -p "Continue? (Y/N): " confirm && [[ "$confirm" == [yY] || "$confirm" == [yY][eE][sS] ]] || exit 1


#if [ "$SETTING" = 1 ];then
#
#  settings used for resfinder criteria selecting nanopore reads for downstream processing
#  settings are conservative so it's expected that every read should have a proper context length for blast analysis
#	CLEANUP=0;
#
#else
#
#  settings used for resfinder criteria selecting nanopore reads for downstream processing
#  settings are opportunistic so it's expected that more reads have a highrer amount of false positives after blast analysis
#	CLEANUP=1;
#
#fi

if [[  "$CONTAMINANT" == H ]] || [[  "$CONTAMINANT" == h ]];then

	DBcont="$dbPATH"/human/;
	shortContaminant=human;
	echo "human sequences will be deleted for downstrean analysis";

elif [[ "$CONTAMINANT" == D ]] || [[ "$CONTAMINANT" == d ]]; then

	DBcont="$dbPATH"/dog/;
	shortContaminant=dog;
	echo "dog sequences will be deleted for downstrean analysis";

elif [[ "$CONTAMINANT" == C ]] || [[ "$CONTAMINANT" == c ]]; then

	DBcont="$dbPATH"/cat/;
	shortContaminant=cat;
	echo "cat sequences will be deleted for downstrean analysis";

elif [[ "$CONTAMINANT" == B ]] || [[ "$CONTAMINANT" == b ]]; then

	DBcont="$dbPATH"/bact/;
	shortContaminant=bact;
	echo "bact sequences will be deleted for downstrean analysis";

elif [[ "$CLEANUP" == 0 ]] ; then

	echo 'no cleanup '
else
 
	echo 'some strange thing happened'
	exit 1

fi


### standard output directories
WORKDIR="$PWD";									                           	# w
CLEANUP="$SETTING";															# x
CONTAMINANT="$CONTAMINANT"													# y
short="$shortContaminant";													# z
RAW_FASTQ="$WORKDIR/RAWREADS/";                                            	# a
RAWSTATS="$WORKDIR/01_rawstats/";                                          	# b
POLISHED="$WORKDIR/02_polished/";                                        	# c
TRIMMEDSTATS="$WORKDIR/03_trimmedstats/";		                           	# d
MAPPING="$WORKDIR/04_mapping/";                                             # m
KRAKEN="$WORKDIR/05_kraken_krona/";											# e
REFERENCES="$WORKDIR/05_references/";										# f
DB="$DBcont";																# h
selectDB="$SETTING2";														# i



BLAST="$WORKDIR/07_blast/";                                                 # h

REPORTING="$WORKDIR/REPORTING/";                                           	# r
TMP="$WORKDIR/TEMP/";                                                      	# s
LOG="$WORKDIR/LOGS/";                                                      	# t





## first steps several scripts with selection steps to reduce data

##  input raw nanopore reads
##  samples have native barcodes so demultiplexing is done bij de gridion


##  structure of the data
./00_structure.sh -w "$WORKDIR" -x "$CLEANUP" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$POLISHED" -d "$TRIMMEDSTATS" -e "$KRAKEN" -f "$REFERENCES" -m "$MAPPING" -r "$REPORTING" -s "$TMP" -t "$LOG"

##  QC of the reads before starting the scripts
#./01_fastqc.sh -w "$WORKDIR" -x "$CLEANUP" -y "$CONTAMINANT" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$POLISHED" -d "$TRIMMEDSTATS" -r "$REPORTING" -s "$TMP" -t "$LOG"


##  trimm reads for artefact
#./02_polish.sh -w "$WORKDIR" -x "$CLEANUP" -y "$CONTAMINANT" -z "$short" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$POLISHED" -d "$TRIMMEDSTATS" -h "$DB" -r "$REPORTING" -s "$TMP" -t "$LOG"


## QC of the reads after trimming the reads for artefact
#./03_fastqctrimmed.sh -w "$WORKDIR" -x "$CLEANUP" -y "$CONTAMINANT" -z "$short" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$POLISHED" -d "$TRIMMEDSTATS" -h "$DB" -r "$REPORTING" -s "$TMP" -t "$LOG"


## align polished reads to viral databases
#./04_mapping_VirCap.sh -w "$WORKDIR" -x "$CLEANUP" -y "$CONTAMINANT" -z "$short" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$POLISHED" -d "$TRIMMEDSTATS" -e "$KRAKEN" -f "$REFERENCES" -m "$MAPPING" -r "$REPORTING" -s "$TMP" -t "$LOG"

## parse output files to construct a top list for mapping and visualization using IGV
#./04p_parse-output.sh -w "$WORKDIR" -x "$CLEANUP" -z "$short" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$POLISHED" -d "$TRIMMEDSTATS" -e "$KRAKEN" -f "$REFERENCES" -m "$MAPPING" -r "$REPORTING" -s "$TMP" -t "$LOG"

## kraken analysis from samples
./05_kraken_krona.sh -w "$WORKDIR" -x "$CLEANUP" -z "$short" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$POLISHED" -d "$TRIMMEDSTATS" -e "$KRAKEN" -f "$REFERENCES" -i "$selectDB" -m "$MAPPING" -r "$REPORTING" -s "$TMP" -t "$LOG"

## resfinder analysis on seperate files, future analysis will be on concat files per sample!
## resfinder_db will be updated by force so the latest will be used for analysis 
#./04_resfinder.sh  -w "$WORKDIR" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$MINLENGTH" -d "$PORECHOP" -e "$FASTA" -f "$RESFINDER" -g "$RESFINDERSELECT" -h "$BLAST" -k "$IDENTITY" -l "$COVERAGE" -m "$CONTEXTLEN" -n "$MINLEN" -r "$REPORTING" -s "$TMP" -t "$LOG" 

## selection of reads from resfinder analysis for downstream balst analysis
#./05_resfinder-read-name.sh  -w "$WORKDIR" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$MINLENGTH" -d "$PORECHOP" -e "$FASTA" -f "$RESFINDER" -g "$RESFINDERSELECT" -h "$BLAST" -k "$IDENTITY" -l "$COVERAGE" -m "$CONTEXTLEN" -n "$MINLEN" -r "$REPORTING" -s "$TMP" -t "$LOG" 

## mask the regions of the fasta files containing the reads with n's so blast analysis will be focused on the context.
#./06_mask-fasta.sh  -w "$WORKDIR" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$MINLENGTH" -d "$PORECHOP" -e "$FASTA" -f "$RESFINDER" -g "$RESFINDERSELECT" -h "$BLAST" -k "$IDENTITY" -l "$COVERAGE" -m "$CONTEXTLEN" -n "$MINLEN" -r "$REPORTING" -s "$TMP" -t "$LOG" 

## blast analysis on the masked reads
#./07_blastn.sh -w "$WORKDIR" -a "$RAW_FASTQ" -b "$RAWSTATS" -c "$MINLENGTH" -d "$PORECHOP" -e "$FASTA" -f "$RESFINDER" -g "$RESFINDERSELECT" -h "$BLAST" -k "$IDENTITY" -l "$COVERAGE" -m "$CONTEXTLEN" -n "$MINLEN" -r "$REPORTING" -s "$TMP" -t "$LOG" 



exit 1
