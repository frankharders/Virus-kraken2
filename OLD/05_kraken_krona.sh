#!/bin/bash

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate NGS-kraken2;

## versions of the used scripts/programs will be documented
PROGVERSIONS="$PWD"/sofware-versions-used.txt;

		ENV=();
		ENV=NGS-kraken2;
		PROG=();
		PROG=$(conda list | grep '^kraken2' | sed 's/  */;/g' | cut -f2 -d';');
		PROG1=();
		PROG1=$(conda list | grep '^krona' | sed 's/  */;/g' | cut -f2 -d';');


			echo -e "script: 05_kraken_krona.sh" >> "$PROGVERSIONS";
			echo -e "conda environment=$ENV" >> "$PROGVERSIONS";			
			echo -e "programm\tversion" >> "$PROGVERSIONS";
			echo -e "kraken2\t$PROG" >> "$PROGVERSIONS";	
			echo -e "krona\t$PROG1" >>  "$PROGVERSIONS";
			echo -e "all scripts used will be documented in this file" >> "$PROGVERSIONS";
			echo -e "\n" >> "$PROGVERSIONS";

while getopts "w:x:y:z:a:b:c:d:e:f:i:m:r:s:t:" opt; do
  case $opt in
     w)
      echo "-w was triggered! $OPTARG"
      WORKDIR="`echo $(readlink -m $OPTARG)`"
      echo $WORKDIR
      ;;
	 x)
      echo "-x was triggered! $OPTARG"
      CLEANUP="$OPTARG"
      echo $CLEANUP
      ;;
	 y)
      echo "-y was triggered! $OPTARG"
      CONTAMINANT="$OPTARG"
      echo $CONTAMINANT
      ;;
	 z)
      echo "-z was triggered! $OPTARG"
      short="$OPTARG"
      echo $short
      ;;	  
     a)
      echo "-a was triggered! $OPTARG"
      RAW_FASTQ="`echo $(readlink -m $OPTARG)`"
      echo $RAW_FASTQ
      ;;
     b)
      echo "-b was triggered! $OPTARG"
      RAWSTATS="`echo $(readlink -m $OPTARG)`"
      echo $RAWSTATS
      ;;
	 c)
      echo "-c was triggered! $OPTARG"
      POLISHED="`echo $(readlink -m $OPTARG)`"
      echo $POLISHED
      ;;
	 d)
      echo "-d was triggered! $OPTARG"
      TRIMMEDSTATS="`echo $(readlink -m $OPTARG)`"
      echo $TRIMMEDSTATS
      ;;
	 e)
      echo "-e was triggered! $OPTARG"
      KRAKEN="`echo $(readlink -m $OPTARG)`"
      echo $KRAKEN
      ;;
	 f)
      echo "-e was triggered! $OPTARG"
      REFERENCES="`echo $(readlink -m $OPTARG)`"
      echo $REFERENCES
      ;;	 
	 i)
      echo "-i was triggered! $OPTARG"
      selectDB="$OPTARG"
      echo $selectDB
      ;;	 
	 m)
      echo "-r was triggered! $OPTARG"
      MAPPING="`echo $(readlink -m $OPTARG)`"
      echo $MAPPING
      ;;
     r)
      echo "-r was triggered! $OPTARG"
      REPORTING="`echo $(readlink -m $OPTARG)`"
      echo $REPORTING
      ;;
     s)
      echo "-s was triggered! $OPTARG"
      TMP="`echo $(readlink -m $OPTARG)`"
      echo $TMP
      ;;
     t)
      echo "-t was triggered! $OPTARG"
      LOG="`echo $(readlink -m $OPTARG)`"
      echo $LOG
      ;;	  
	\?)
      echo "$OPTARG" >&2
      ;;

  esac
done



cnt=$(cat samples.txt | wc -l);



##  variables
CPU=24;
CL=0;
MRC=10;
KL=150;
threshold="$MRC";
LEVEL=S;

HEAD="$TMP"/header.temp.bracken.tab;
FILES="$TMP"/concat.bracken.outputs.temp.tab;

# remove old concatenated file
rm "$REPORTING"/all.samples.compare.bracken;
rm "$FILES";
rm "$HEAD";
rm "$REPORTING"/all.samples.c."$CL".db."$DB".level."$LEVEL".compare.final.bracken;



count0=1;
countS=$(cat samples.txt | wc -l);

while [ "$count0" -le "$countS" ]; do 

	SAMPLE=$(cat samples.txt | awk 'NR=='"$count0");



if [ "$CLEANUP" == 1 ];then
	R1="$POLISHED"/"$SAMPLE"_R1.QTR.adapter.correct.fq.gz;
	R2="$POLISHED"/"$SAMPLE"_R2.QTR.adapter.correct.fq.gz;
else
	R1="$POLISHED"/"$SAMPLE"_R1.QTR."$short".clean.adapter.correct.fq.gz;
	R2="$POLISHED"/"$SAMPLE"_R2.QTR."$short".clean.adapter.correct.fq.gz;
fi



## select DB for analysis & plots

if [ "$selectDB" = 0 ] ; then

	DBpath='/mnt/lely_DB/kraken2/k2_standard_20221209';
	DB='k2_standard';

elif [ "$selectDB" = 1 ];then

	DBpath='/mnt/lely_DB/kraken2/k2_viral_20221209';
	DB='k2_virus';

elif [ "$selectDB" = 2 ]; then

	DBpath='/mnt/lely_DB/kraken2/kraken_DB_08092022_pluspfp';
	DB='k2_pluspfp';
else

exit 1
	
fi

##  output files
KRAKENout1="$KRAKEN"/"$SAMPLE".c."$CL".db."$DB".tsv;
KRAKENout2="$KRAKEN"/"$SAMPLE".c."$CL".db."$DB".report.tsv;
KRONA1="$KRAKEN"/"$SAMPLE".c."$CL".db."$DB".krona;
KRONA2="$KRAKEN"/"$SAMPLE".c."$CL".db."$DB".html;
KRAKENout3="$KRAKEN"/"$SAMPLE".1.c."$CL".db."$DB".classified.fasta;
KRAKENout4="$KRAKEN"/"$SAMPLE".2.c."$CL".db."$DB".classified.fasta;
BRACKENout="$KRAKEN"/"$SAMPLE".c."$CL".db."$DB".kmers.$KL.level."$LEVEL".abundance.bracken.tab;



##  kraken2 analysis
kraken2 --db "$DBpath" --use-names --report-minimizer-data --report "$KRAKENout2" --threads "$CPU" --confidence "$CL" --output "$KRAKENout1" --paired --classified-out "$KRAKEN"/"$SAMPLE".db."$DB".cseqs#.fq "$R1" "$R2";



##  kraken 2 krona plots
/home/harde004/GIT/KrakenTools/kreport2krona.py -r "$KRAKENout2" -o "$KRONA1";

## krona sunburst plots in html format
ktImportText "$KRONA1" -o "$KRONA2"; 


#NAMES=/mnt/lely_DB/kraken2/taxonomy/names.dmp;
#NODES=/mnt/lely_DB/kraken2/taxonomy/nodes.dmp;

# problems using the taxonomy in S is used and the samples conatiing S1 and S2 taxid's
#stringmeup --names $NAMES --nodes $NODES 0.1 $KRAKENout1 --output_classifications $KRAKEN/$SAMPLE.class.new.tab;

d MY_DB -i INPUT -o OUTPUT -w OUTREPORT -r READ_LEN -l LEVEL -t THRESHOLD

~/GIT/Bracken/bracken -d "$DBpath" -i "$KRAKENout2" -o default -w "$BRACKENout" -r 150 -l "$LEVEL" -t "$MRC";

#python ~/GIT/Bracken/src/est_abundance.py -i "$KRAKENout2" -k "$DBpath"/'database'$KL'mers.kmer_distrib' -l S -t "$threshold" -o "$BRACKENout";


## extract fasta from the list above for mapping purposes
#while read i;do extract_kraken_reads.py -k ./05_kraken_krona/1.c.0.db.virus.tsv -s ./02_polished/1_R1.QTR.adapter.correct.fq.gz -s2 ./02_polished/1_R2.QTR.adapter.correct.fq.gz  -t $i  -o #./02_polished/1_R1.kraken.$i.fastq -o2 ./02_polished/1_R2.kraken.$i.fastq -r ./05_kraken_krona/1.c.0.db.virus.report.tsv --fastq-output ; done < ./05_kraken_krona/test.names.output5


count1=2; # because the file contains a 1 line header
countT=$(cat "$BRACKENout" | wc -l);

cat "$BRACKENout" | head -n1 > "$HEAD";

while [ "$count1" -le "$countT" ];do

LINE=$(cat "$BRACKENout" | awk 'NR=='"$countT");

#echo $LINE;

	echo -e "$SAMPLE\t$LINE" >> "$FILES";

count1=$((count1+1));
done < "$BRACKENout"


	let cnt--;	
	echo -e "$cnt samples to go!";
	echo "NEXT";

count0=$((count0+1));

done

cat "$HEAD" "$FILES" > "$REPORTING"/all.samples.c."$CL".db."$DB".level."$LEVEL".compare.final.bracken;



## eerste 2 lijnen van het report in col 2 staan opgesomd het aantal reads in de analyse
#
#TOTreadcnt=$(cat "$KRAKENout2" | head -n2 | cut -f2 -d$'\t' |paste -sd+ | bc);
#
#cat "$KRAKENout2" | grep -nP '\tD\t' | cut -f5,6 -d$'\t' > "$TMP"/"$SAMPLE".Dlist;
#
#while read LINE;do 
#col2=$(echo "$LINE" | cut -f2 -d$'\t');
#col3=$(echo "$LINE" | cut -f3);
#col4=$(echo "$LINE" | cut -f4);
#col5=$(echo "$LINE" | cut -f5);
#col6=$(echo "$LINE" | cut -f6);
#
#
#if [ "$col2" -gt 10 ] && [ "$col2" -eq "$col3" ];then 
#
#
#echo -e "$col2\t$col3\t$col4\t$col5\t$col6" >> test.names.output;
#
#fi     
#
#done < 9.c.0.db.virus.report.tsv









##### old! #####
#krakenDB='/mnt/lely_DB/KRAKEN2/standard-prot-fungi';
#CONFIDENCE=0.95;
#NODES=48;
#TAX='/mnt/lely_DB/Taxonomy/';
#
#cnt=$(cat samples.txt | wc -l);
#
#
#count0=1;
#countS=$(cat samples.txt | wc -l);
#
#while [ "$count0" -le "$countS" ]; do 
#
#	SAMPLE=$(cat samples.txt | awk 'NR=='"$count0");
#
#LOG1="$LOG"/"$SAMPLE".general.log;
#
#echo "fastQC analysis on trimmed sequence files";
#
#if [ "$CLEANUP" == 1 ];then
#	R1="$POLISHED"/"$SAMPLE"_R1.QTR.adapter.correct.fq.gz;
#	R2="$POLISHED"/"$SAMPLE"_R2.QTR.adapter.correct.fq.gz;
#
#		OUT="$KRAKEN"/"$SAMPLE".paired.tsv;
#		OUTreport="$KRAKEN"/"$SAMPLE".paired.report.tsv;
#		OUThtml="$KRAKEN"/"$SAMPLE".paired.html;
#
# kraken analysis
#			kraken2 --db "$krakenDB" --paired "$R1" "$R2" --confidence "$CONFIDENCE" --output "$OUT" --report "$OUTreport" --threads "$NODES"  >> "$LOG1" 2>&1;
#
# create html for visualisation
#				cut -f2,3 "$OUTreport" | ktImportTaxonomy -k -o "$OUThtml" -tax "$DBkronanew" -  >> "$LOG1" 2>&1;
#
#
#else
#	R1="$POLISHED"/"$SAMPLE"_R1.QTR."$short".clean.adapter.correct.fq.gz;
#	R2="$POLISHED"/"$SAMPLE"_R2.QTR."$short".clean.adapter.correct.fq.gz;
#
#		OUT="$KRAKEN"/"$SAMPLE"."$short".clean.paired.tsv;
#		OUTreport="$KRAKEN"/"$SAMPLE"."$short".clean.paired.report.tsv;
#		OUThtml="$KRAKEN"/"$SAMPLE"."$short".clean.paired.html;
#
# kraken analysis
#			kraken2 --db "$krakenDB" --paired "$R1" "$R2" --confidence "$CONFIDENCE" --output "$OUT" --report "$OUTreport" --threads "$NODES"  >> "$LOG1" 2>&1;
#
# create html for visualisation
#				cut -f2,3 "$OUTreport" | ktImportTaxonomy -k -o "$OUThtml" -tax "$DBkronanew" -  >> "$LOG1" 2>&1;
#
#fi
#
#	let cnt--;	
#	echo -e "$cnt samples to go!";
#	echo "NEXT";
#
#count0=$((count0+1));
#
#done

exit 1


##### kraken2
#
#
#Need to specify input filenames!
#Usage: kraken2 [options] <filename(s)>
#
#Options:
#  --db NAME               Name for Kraken 2 DB
#                          (default: none)
#  --threads NUM           Number of threads (default: 1)
#  --quick                 Quick operation (use first hit or hits)
#  --unclassified-out FILENAME
#                          Print unclassified sequences to filename
#  --classified-out FILENAME
#                          Print classified sequences to filename
#  --output FILENAME       Print output to filename (default: stdout); "-" will
#                          suppress normal output
#  --confidence FLOAT      Confidence score threshold (default: 0.0); must be
#                          in [0, 1].
#  --minimum-base-quality NUM
#                          Minimum base quality used in classification (def: 0,
#                          only effective with FASTQ input).
#  --report FILENAME       Print a report with aggregrate counts/clade to file
#  --use-mpa-style         With --report, format report output like Kraken 1's
#                          kraken-mpa-report
#  --report-zero-counts    With --report, report counts for ALL taxa, even if
#                          counts are zero
#  --report-minimizer-data With --report, report minimizer and distinct minimizer
#                          count information in addition to normal Kraken report
#  --memory-mapping        Avoids loading database into RAM
#  --paired                The filenames provided have paired-end reads
#  --use-names             Print scientific names instead of just taxids
#  --gzip-compressed       Input files are compressed with gzip
#  --bzip2-compressed      Input files are compressed with bzip2
#  --minimum-hit-groups NUM
#                          Minimum number of hit groups (overlapping k-mers
#                          sharing the same minimizer) needed to make a call
#                          (default: 2)
#  --help                  Print this message
#  --version               Print version information
#
#If none of the *-compressed flags are specified, and the filename provided
#is a regular file, automatic format detection is attempted.
#
#
#####

##### bracken
#Usage: bracken -d MY_DB -i INPUT -o OUTPUT -w OUTREPORT -r READ_LEN -l LEVEL -t THRESHOLD
#  MY_DB          location of Kraken database
#  INPUT          Kraken REPORT file to use for abundance estimation
#  OUTPUT         file name for Bracken default output
#  OUTREPORT      New Kraken REPORT output file with Bracken read estimates
#  READ_LEN       read length to get all classifications for (default: 100)
#  LEVEL          level to estimate abundance at [options: D,P,C,O,F,G,S,S1,etc] (default: S)
#  THRESHOLD      number of reads required PRIOR to abundance estimation to perform reestimation (default: 0)
#####

