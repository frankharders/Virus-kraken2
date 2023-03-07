#!/bin/bash

cnt=$(cat samples.txt | wc -l);

while getopts "i:o:r:w:" opt; do
  case $opt in
    i)
      echo "-i was triggered! $OPTARG"
      SEQDIR="`echo $(readlink -m $OPTARG)`"
      echo $SEQDIR
      ;;
    o)
      echo "-o was triggered! $OPTARG"
      OUTDIR="`echo $(readlink -m $OPTARG)`"
	  echo $OUTDIR
      ;;
	r)
	  echo "-l was triggered! $OPTARG"
	  REFDIR="`echo $(readlink -m $OPTARG)`"
	  echo $REFDIR
	  ;;
	w)
	  echo "-p was triggered! $OPTARG"
	  WORKDIR="`echo $(readlink -m $OPTARG)`"
	  echo $WORKDIR
	  ;;
    \?)
      echo "-i for the folder containing fastq files, -o for output folder of polished data: -$OPTARG" >&2
      ;;
  esac
done

if [ "x" == "x$SEQDIR" ] || [ "x" == "x$OUTDIR" ] || [ "x" == "x$REFDIR" ] || [ "x" == "x$WORKDIR" ] ; then
    echo "-i $SEQDIR -o $OUTDIR -r $REF -w $WORKDIR"
    echo "-i, -o -r -w [options] are required"
    exit 1
fi

# build reference index files

reference='NC_018464.1';
REF=$reference'.fna';

NODES=16;


SAMPLESreadcnt='./REPORTING/mappedreads.'$reference'.tab';

rm $SAMPLESreadcnt;

echo -e "SampleName\tMappedReads" > $SAMPLESreadcnt;# header


while read SAMPLE; do 
mkdir -p $OUTDIR/$SAMPLE/;

MAPDIR=$OUTDIR/$SAMPLE/;
MapPEin1=$SEQDIR/$SAMPLE'_R1.fastq';
MapPEin2=$SEQDIR/$SAMPLE'_R2.fastq';
#MapSEin=$SEQDIR/$SAMPLE'_SE.trimmed.fastq';
SAMout=$OUTDIR/$SAMPLE'_mapped.'$reference'.sam';
LOG1=./LOGS/$SAMPLE'.bowtie2.mapped.'$reference'.log';
BAMout=$OUTDIR/$SAMPLE'_mapped.'$reference'.bam';
BAMsorted=$MAPDIR/$SAMPLE'_mapped.sorted.'$reference'.bam';
BAMbai=$MAPDIR/$SAMPLE'_mapped.sorted.'$reference'.bai';
BOWTIEMETRICS=$MAPDIR/$SAMPLE'_'$REFSHORT'_metrics.mapping.log';
FLAGSTATout=$MAPDIR/$SAMPLE'_'$REFSHORT'.flagstat.tab';
COVERAGEout=$MAPDIR/$SAMPLE'_'$REFSHORT'.coverage.tab';
IDXSTATSout=$MAPDIR/$SAMPLE'_'$REFSHORT'.idxstats.tab';



echo -e "ReferenceName\tGenomeLength\tMappedReadSegments\tUnmappedReadSegments" > $COVERAGEout;


bowtie2 --local --very-sensitive-local --no-unal --phred33 -p 8 -x $REFDIR/$reference -1 $MapPEin1 -2 $MapPEin2 -S $SAMout > $LOG1 2>&1;

		samtools view -@ $NODES -S $SAMout -b -o $BAMout;
		samtools sort -@ $NODES  -o $BAMsorted $BAMout;
		samtools index $BAMsorted $BAMbai;
		#samtools idxstats $BAMsorted >> $IDXSTATSout;
		#samtools flagstat $BAMsorted -o tsv > $FLAGSTATout;
		#samtools coverage $BAMsorted -o $COVERAGEout;

		READcnt=$(samtools view -q 15 $SAMout | cut -f1 | sort | uniq | wc -l ); # count reads directly from sam file
		READcnt2=$(samtools view -q 15 $SAMout | cut -f1 | sort | wc -l ); # count reads directly from sam file

		echo -e "$SAMPLE\t$READcnt\t$READcnt2";
		echo -e "$SAMPLE\t$READcnt" >> $SAMPLESreadcnt;
		

#rm $SAMout $BAMout;



	let cnt--;
	echo -e "$cnt samples to go!";
	echo "NEXT";

		done < samples.txt


exit 1

