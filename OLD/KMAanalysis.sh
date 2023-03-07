#!/bin/bash

KMAout='./KMAoutput/';

mkdir -p $KMAout;

FILEdir='./02_polished/';


count0=1;
countS=$(cat samples.txt | wc -l);


while [ $count0 -le $countS ];do

	SAMPLE=$(cat samples.txt | awk 'NR=='$count0 );



R1=$FILEdir/$SAMPLE/$SAMPLE'_R1.cleaned.QTR.adapter.correct.fq.gz';
R2=$FILEdir/$SAMPLE/$SAMPLE'_R2.cleaned.QTR.adapter.correct.fq.gz';
#REF='VirCapSeq_targets';
REF='viral_20200519';
FILEout=$KMAout/$SAMPLE'.kma.'$REF'.tab';
DBin='/mnt/lely_DB/minikraken2/viral_20200519/viral_20200519_kma/viral_20200519_kma';
#DBin='./references/VirCap/kma_index_vircap/VirCapSeq_targets';

SAMout=$KMAout/$SAMPLE'.kma.'$REF'.sam';
MP=20;
threads=8;

echo $SAMPLE;

kma -t $threads -ipe $R1 $R2 -o $FILEout -t_db $DBin  -nc -mp $MP -fpm p  -ConClave 2 -cge -dense -ef -1t1 -ca -boot;# -sam;

#kma -t $threads -i $R2 -o $FILEout -t_db $DBin  -nc -mp $MP -fpm p  -ConClave 2 -cge -dense -ef -1t1 -ca -boot;# -sam;


count0=$((count0+1));

done

exit 1



#KMA manual
# KMA-1.2.3 mapps raw reads to a template database.
# Options are:		Desc:				Default:	Requirements:
#
#	-o		Output file			None		REQUIRED
#	-t_db		Template DB			None		REQUIRED
#	-i		Input file name(s)		STDIN
#	-ipe		Input paired end file name(s)
#	-int		Input interleaved file name(s)
#	-k		Kmersize			DB defined
#	-e		evalue				0.05
#	-ConClave	ConClave version		1
#	-mem_mode	Use kmers to choose best
#			template, and save memory	False
#	-ex_mode	Searh kmers exhaustively	False
#	-ef		Print additional features	False
#	-vcf		Make vcf file, 2 to apply FT	False/0
#	-sam		Output sam, 4 to only output
#			mapped reads			False/0
#	-nc		No consensus file		False
#	-nf		No frag file		False
#	-deCon		Remove contamination		False
#	-dense		Do not allow insertions
#			in assembly			False
#	-ref_fsa	Consensus sequnce will
#			have "n" instead of gaps	False
#	-matrix		Print assembly matrix		False
#	-a		Print all best mappings		False
#	-mp		Minimum phred score		20
#	-5p		Cut a constant number of
#			nucleotides from the 5 prime.	0
#	-Sparse		Only count kmers		False
#	-Mt1		Map only to "num" template.	0 / False
#	-ID		Minimum ID			1.0%
#	-ss		Sparse sorting (q,c,d)		q
#	-pm		Pairing method (p,u,f)		u
#	-fpm		Fine Pairing method (p,u,f)	u
#	-apm		Sets both pm and fpm		u
#	-shm		Use shared DB made by kma_shm	0 (lvl)
#	-1t1		Force end to end mapping	False
#	-ck		Count kmers instead of
#			pseudo alignment		False
#	-ca		Make circular alignments	False
#	-boot		Bootstrap sequence		False
#	-bc		Base calls should be
#			significantly overrepresented.	[True]
#	-bc90		Base calls should be both
#			significantly overrepresented,
#			and have 90% agreement.		False
#	-bcNano		Call bases at suspicious
#			deletions, made for nanopore.	False
#	-bcd		Minimum depth at base		1
#	-bcg		Maintain insignificant gaps
#	-and		Both mrs and p_value thresholds
#			has to reached to in order to
#			report a template hit.		or
#	-mq		Minimum mapping quality		0
#	-mrs		Minimum alignment score,
#			normalized to alignment length	0.50
#	-reward		Score for match			1
#	-penalty	Penalty for mismatch		-2
#	-gapopen	Penalty for gap opening		-3
#	-gapextend	Penalty for gap extension	-1
#	-per		Reward for pairing reads	7
#	-cge		Set CGE penalties and rewards	False
#	-t		Number of threads		1
#	-v		Version
#	-h		Shows this help message
#
