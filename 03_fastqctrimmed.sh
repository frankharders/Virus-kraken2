#!/bin/bash

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate NGS-QC;

## versions of the used scripts/programs will be documented
PROGVERSIONS="$PWD"/sofware-versions-used.txt;

		ENV=();
		ENV=NGS-QC;
		PROG=();
		PROG=$(conda list | grep '^fastqc' | sed 's/  */;/g' | cut -f2 -d';');

			echo -e "script: 03_fastqctrimmed.sh" >> "$PROGVERSIONS";
			echo -e "conda environment=$ENV" >> "$PROGVERSIONS";			
			echo -e "programm\tversion" >> "$PROGVERSIONS";
			echo -e "fastqc\t$PROG" >> "$PROGVERSIONS";	
			echo -e "all scripts used will be documented in this file" >> "$PROGVERSIONS";
			echo -e "\n" >> "$PROGVERSIONS";


while getopts "w:x:y:z:a:b:c:d:f:h:m:r:s:t:" opt; do
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


count0=1;
countS=$(cat samples.txt | wc -l);

while [ "$count0" -le "$countS" ]; do 

	SAMPLE=$(cat samples.txt | awk 'NR=='"$count0");

LOG1="$LOG"/"$SAMPLE".general.log;

echo "fastQC analysis on trimmed sequence files";

if [ "$CLEANUP" == 1 ];then
	R1="$POLISHED"/"$SAMPLE"_R1.QTR.adapter.correct.fq.gz;
	R2="$POLISHED"/"$SAMPLE"_R2.QTR.adapter.correct.fq.gz;

		fastqc -t 8 -o $TRIMMEDSTATS $R1 $R2 >> "$LOG1" 2>&1;

else
	R1="$POLISHED"/"$SAMPLE"_R1.QTR."$short".clean.adapter.correct.fq.gz;
	R2="$POLISHED"/"$SAMPLE"_R2.QTR."$short".clean.adapter.correct.fq.gz;

		fastqc -t 8 -o $TRIMMEDSTATS $R1 $R2 >> "$LOG1" 2>&1;

fi

	let cnt--;	
	echo -e "$cnt samples to go!";
	echo "NEXT";

count0=$((count0+1));

done

rm "$TRIMMEDSTATS"/*.zip;


echo "fastqc plots are generated for all polished reads for all samples";
echo -e "output files for downstream processing can be found in the directory $TRIMMEDSTATS";



	
exit 1

