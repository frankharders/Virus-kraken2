#!/bin/bash

# create directories
cnt=$(cat samples.txt | wc -l);



count0=1;
countf=$(cat samples.txt | wc -l);

while [ $count0 -le $countf ];do

	SAMPLE=$(cat samples.txt | awk 'NR=='$count0 );

echo $SAMPLE;

KMAres=./KMAoutput/$SAMPLE'.kma.tab.res';
KMAresout=./KMAoutput/$SAMPLE'.kma.tab.res.sample';



countL=1;
countLF=$(cat $KMAres | wc -l);

while [ $countL -le $countLF ];do

LINE=$(cat $KMAres | awk 'NR=='$countL);

echo -e "$SAMPLE\t$LINE" >> $KMAresout;



countL=$((countL+1));

done 

KMAmapstat=./KMAoutput/$SAMPLE'.kma.tab.mapstat';
KMAmapstatout=./KMAoutput/$SAMPLE'.kma.tab.mapstat.sample';


countM=1;
countMF=$(cat $KMAmapstat | wc -l);

while [ $countM -le $countMF ];do

LINES=$(cat $KMAmapstat | awk 'NR=='$countM);

echo -e "$SAMPLE\t$LINES" >> $KMAmapstatout;


countM=$((countM+1));

done 

count0=$((count0+1));

done 

cat ./KMAoutput/*.kma.tab.res.sample > ./KMAoutput/allsamples.kma.res.tab;
cat ./KMAoutput/*.kma.tab.mapstat.sample > ./KMAoutput/allsamples.kma.mapstat.tab;



exit 1






















