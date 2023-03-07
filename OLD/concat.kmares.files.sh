#!/bin/bash

count0=1;
countS=$(cat samples.txt | wc -l);

while [ "$count0" -le "$countS" ]; do 

	SAMPLE=$(cat samples.txt | awk 'NR=='"$count0");


MAPP="$PWD"/04_mapping;
MAPPout1="$PWD"/04_mapping/all.samples.vircap.tab.res;
MAPPout2="$PWD"/04_mapping/all.samples.ncbidatabase.tab.res;


while read LINE;do

echo -e "$SAMPLE\t$LINE" >> "$MAPPout1";




done < "$MAPP"/"$SAMPLE".vircap.tab.res




while read LINE;do

echo -e "$SAMPLE\t$LINE" >> "$MAPPout2";




done < "$MAPP"/"$SAMPLE".viral_20200519.tab.res



MAPP="$PWD"/04_mapping;
MAPPout3="$PWD"/04_mapping/all.samples.vircap.tab.mapstat;
MAPPout4="$PWD"/04_mapping/all.samples.ncbidatabase.tab.mapstat;



while read LINE;do

echo -e "$SAMPLE\t$LINE" >> "$MAPPout3";




done < "$MAPP"/"$SAMPLE".vircap.tab.mapstat




while read LINE;do

echo -e "$SAMPLE\t$LINE" >> "$MAPPout4";




done < "$MAPP"/"$SAMPLE".viral_20200519.tab.mapstat










count0=$((count0+1));
done

exit 1
