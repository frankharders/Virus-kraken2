#!/bin/bash

mkdir -p ./00_cleaned;
mkdir -p ./LOGS;

SAMPLEFILE='samples.txt';

count0=1;
countS=$(cat $SAMPLEFILE | wc -l);

dbPATH='/mnt/lely_DB/Old_Shared_DATABASES/CONTAMINANTSCREENING/human/';


while [ $count0 -le $countS ];do

SAMPLE=$(cat $SAMPLEFILE | awk 'NR=='$count0);

FILEin1=./RAWREADS/$SAMPLE'_R1.fastq';
FILEin2=./RAWREADS/$SAMPLE'_R2.fastq';
FILEout1=./00_cleaned/$SAMPLE'_R1.cleaned.fq.gz';
FILEout2=./00_cleaned/$SAMPLE'_R2.cleaned.fq.gz';



bbmap.sh -Xmx64g minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=$dbPATH qtrim=rl trimq=10 untrim in1=$FILEin1 in2=$FILEin2 outu1=$FILEout1 outu2=$FILEout2 outm=human.fq ow=t > ./LOGS/$SAMPLE'.cleaned.log' 2>&1;


count0=$((count0+1));

done

exit 1
