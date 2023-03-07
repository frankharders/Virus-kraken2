#!/bin/bash


OUT_FOLDER="$PWD";
SAMPLEFILE='samples.txt';

for REF in `ls $OUT_FOLDER/references/resfinder_old/*.fsa`;do

	REFSHORT=`echo $REF | cut -f1 -d'.' | cut -f11 -d'/'`;
	echo -e "$REFSHORT";

while read SAMPLE;do

	RESOUT="$OUT_FOLDER/03_mapping/mapping_$REFSHORT";
	outdir=$RESOUT/$SAMPLE/;
#	RPKM=$outdir/$SAMPLE'.'$REFSHORT'.RPKM.tab';
#	REFSTATS=$outdir/$SAMPLE'.'$REFSHORT'.refstats.tab';
	SCAFSTATS=$outdir/$SAMPLE'.'$REFSHORT'.scafstats.tab';

	parse1=$outdir/$SAMPLE'.scafstats.'$REFSHORT'.head.txt';
	parse2=$outdir/$SAMPLE'.scafstats.'$REFSHORT'.headless.txt';
	parse3=$outdir/$SAMPLE'.scafstats.'$REFSHORT'.gene.headless.txt';
	parse4=$outdir/all_samples.scafstats.headless.$REFSHORT.tab;

rm $parse3;
#rm $parse4;



	head=$(cat $SCAFSTATS | head -n1);
	echo -e "SAMPLE\tREFSHORT\tGENE\t$head" > $parse1;
	cat $SCAFSTATS | sed 1d > $parse2;


while read LINE;do

		seq=$(echo "$LINE" | cut -f1 -d$'\t');
		gene=$(echo $seq | cut -f1 -d'_' );
		#genesupershort=$(echo "$gene" | cut -f1 -d'_' | cut -c1-3) ;
		#geneshort=$(echo "$gene" | cut -f1 -d'_');
		#class=$(cat $ABclass | grep ""^$genesupershort"" | cut -f2 -d':' | cut -f1 -d' ' | head -n1 | cut -f1 -d',');
	
		#echo -e "$SAMPLE\t$LINE\t$gene\t$geneshort\t$class" ;#>> $parse3;
		#echo -e "$SAMPLE\t$LINE\t$gene\t$geneshort\t$class" ;#>> $parse4;

echo -e "$SAMPLE\t$REFSHORT\t$gene\t$LINE" >> $parse3;



done < $parse2

cat $parse1 $parse3 > $outdir/$SAMPLE'.'$REFSHORT'.final.scafstats.tab';


done < $SAMPLEFILE # end loop sample list

rm all_samples.final.tab;

for k in `ls $OUT_FOLDER/03_mapping/*/*/*.final.scafstats.tab `;do

#echo $k;

cat $k >> all_samples.final.tab;


done

done # end loop ref files --> moet of fsa staan ipv fasta --> testing met fasta file!

cat all_samples.final.tab | sort -Vu > all_samples.sorted.final.tab; 



exit 1


