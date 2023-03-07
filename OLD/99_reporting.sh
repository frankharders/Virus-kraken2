#!/bin/bash


while getopts "w:n:l:r:" opt; do
  case $opt in
       w)
      echo "-w was triggered! $OPTARG"
      WORKDIR="`echo $(readlink -m $OPTARG)`"
      echo $WORKDIR
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
      SHOVILL="`echo $(readlink -m $OPTARG)`"
      echo $SHOVILL
      ;;
	 f)
      echo "-f was triggered! $OPTARG"
      QUAST="`echo $(readlink -m $OPTARG)`"
      echo $QUAST
      ;;
	 g)
      echo "-g was triggered! $OPTARG"
      QUASTparse="`echo $(readlink -m $OPTARG)`"
      echo $QUASTparse
      ;;
     h)
      echo "-h was triggered! $OPTARG"
      MLST="`echo $(readlink -m $OPTARG)`"
      echo $MLST
      ;;
     i)
      echo "-i was triggered! $OPTARG"
      MLSTparse="`echo $(readlink -m $OPTARG)`"
      echo $MLSTparse
      ;;	  
	 j)
      echo "-j was triggered! $OPTARG"
      ABRICATE="`echo $(readlink -m $OPTARG)`"
      echo $ABRICATE
      ;;
     k)
      echo "-k was triggered! $OPTARG"
      ABRICATEparse="`echo $(readlink -m $OPTARG)`"
      echo $ABRICATEparse
      ;;	  
	 l)
      echo "-l was triggered! $OPTARG"
      TMP="`echo $(readlink -m $OPTARG)`"
      echo $TMP
      ;;	  
	 n)
      echo "-n was triggered! $OPTARG"
      LOG="`echo $(readlink -m $OPTARG)`"
      echo $LOG
      ;;
     r)
      echo "-r was triggered! $OPTARG"
      REPORTING="`echo $(readlink -m $OPTARG)`"
      echo $REPORTING
      ;;
	\?)
      echo "-i for the folder containing fastq files, -o for output folder of polished data: -$OPTARG" >&2
      ;;
  esac
done


if [ "x" == "x$WORKDIR" ] || [ "x" == "x$LOG" ] || [ "x" == "x$TMP" ] || [ "x" == "x$REPORTING" ]; then
    echo "-w $WORKDIR -n $LOG -l $TMP -r $REPORTING"
    echo "-w, -n, -l, -r  [options] are required"
    exit 1
fi



METRICSout=$REPORTING/'all.samples.general.metrics.tab';
HEADER=$TMP/metrics.header.txt;
METRICStmp=$TMP/metrics.txt;

rm $METRICStmp $HEADER $METRICSout;

echo -e "Sample\testGenomeSize\tcontigs\tInserSize\tN50\tGC\treads\tmappedReads\tavgCoverage" > $HEADER;
echo -e "Sample\testGenomeSize\tcontigs\tInserSize\tN50\tGC\treads\tmappedReads\tavgCoverage";

while read SAMPLE;do

cp $WORKDIR/01_rawstats/$SAMPLE/*.html $REPORTING/$SAMPLE/;
cp $WORKDIR/03_trimmedstats/$SAMPLE/*.html $REPORTING/$SAMPLE/;
cp $WORKDIR/04_shovill/$SAMPLE/$SAMPLE'_contigs.fa' $REPORTING/$SAMPLE/;
cp $WORKDIR/05_quast_analysis/$SAMPLE/$SAMPLE'.report.pdf' $REPORTING/$SAMPLE/;
cp $WORKDIR/05_quast_analysis/$SAMPLE/$SAMPLE'.report.tsv' $REPORTING/$SAMPLE/;
cp $WORKDIR/06_mlst/$SAMPLE'.mlst.csv' $REPORTING/$SAMPLE/;


InsertMetrics=$WORKDIR/02_polished/$SAMPLE/$SAMPLE'.ihist.txt';
InsertSize=$(cat $InsertMetrics | grep "#Median"| cut -f2);

#tsv file parsing
QUASTin=$WORKDIR/05_quast_analysis/$SAMPLE/transposed_report.tsv;
QUASTparse=$(cat $QUASTin | sed '1d');
N50=$(cat $QUASTin | grep "^$SAMPLE"'_contigs' | cut -f18 );
estGenomeSize=$(cat $QUASTin | grep "^$SAMPLE"'_contigs'| cut -f16 );
contigs=$(cat $QUASTin | grep "^$SAMPLE"'_contigs'| cut -f2);
GC=$(cat $QUASTin | grep "^$SAMPLE"'_contigs'| cut -f17 );
reads=$(cat $QUASTin | grep "^$SAMPLE"'_contigs'| cut -f22 );
mappedReads=$(cat $QUASTin | grep "^$SAMPLE"'_contigs'| cut -f25 );
avgCoverage=$(cat $QUASTin | grep "^$SAMPLE"'_contigs'| cut -f27 );


echo -e "$SAMPLE\t$estGenomeSize\t$contigs\t$InsertSize\t$N50\t$GC\t$reads\t$mappedReads\t$avgCoverage" >> $METRICStmp;
echo -e "$SAMPLE\t$estGenomeSize\t$contigs\t$InsertSize\t$N50\t$GC\t$reads\t$mappedReads\t$avgCoverage";


done < samples.txt

cat $HEADER $METRICStmp > $METRICSout;




exit 1
