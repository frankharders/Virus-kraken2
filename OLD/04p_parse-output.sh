#!/bin/bash

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate NGS-bbmap;

## versions of the used scripts/programs will be documented
PROGVERSIONS="$PWD"/sofware-versions-used.txt;

		ENV=();
		ENV=NGS-bbmap;
		PROG=();
		PROG=$(conda list | grep '^bbmap' | sed 's/  */;/g' | cut -f2 -d';');
		PROG1=();
		PROG1=$(conda list | grep '^samtools' | sed 's/  */;/g' | cut -f2 -d';');


			echo -e "script: 04p_parse-output.sh" >> "$PROGVERSIONS";
			echo -e "conda environment=$ENV" >> "$PROGVERSIONS";			
			echo -e "programm\tversion" >> "$PROGVERSIONS";
			echo -e "bbmap-suite\t$PROG" >> "$PROGVERSIONS";
			echo -e "samtools\t$PROG1" >> "$PROGVERSIONS";

while getopts "w:x:z:a:b:c:d:e:f:m:r:s:t:" opt; do
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
      echo "-d was triggered! $OPTARG"
      REFERENCES="`echo $(readlink -m $OPTARG)`"
      echo $REFERENCES
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


SELECT=25;


# create directories
cnt=$(cat samples.txt | wc -l);

count0=1;
countS=$(cat samples.txt | wc -l);

while [ "$count0" -le "$countS" ]; do 

	SAMPLE=$(cat samples.txt | awk 'NR=='"$count0");


# output files
	RPKM1="$MAPPING"/"$SAMPLE"."$DBshort1".RPKM.tab;
	COVSTATS1="$MAPPING"/"$SAMPLE"."$DBshort1".covstats.tab;
	OUTM1="$TMP"/"$SAMPLE"."$DBshort1".sam;
	OUTMb1="$TMP"/"$SAMPLE"."$DBshort1".bam;

	RPKM2="$MAPPING"/"$SAMPLE"."$DBshort2".RPKM.tab;
	COVSTATS2="$MAPPING"/"$SAMPLE"."$DBshort2".covstats.tab;
	OUTM2="$TMP"/"$SAMPLE"."$DBshort2".sam;
	OUTMb2="$TMP"/"$SAMPLE"."$DBshort2".bam;

if [ "$CLEANUP" == 1 ];then

# variables
TAGS='mdtag=t nhtag=t xmtag=t amtag=t nmtag=t xstag=t';
LOG1="$LOG"/"$SAMPLE".general.log;
DBshort1='viral_20200519';
DBshort2='vircap';
DB1='/mnt/lely_DB/minikraken2/viral_20200519/library/viral/library.fna';
DB2='/mnt/lely_DB/Old_Shared_DATABASES/VIRUSES/VirCap/VirCapSeq_targets.fa.gz';

# select a dedicated number of positives from the bbmap and kma output files from the different databases  

# all input and output kma files
KMALIST1="$MAPPING"/"$SAMPLE".kma."$DBshort1".names;
KMALIST2="$MAPPING"/"$SAMPLE".kma."$DBshort2".names;
KMAin1=/mnt/lely_DB/minikraken2/viral_20200519/library/viral/library.fna;
KMAin2=/mnt/lely_DB/Old_Shared_DATABASES/VIRUSES/VirCap/VirCapSeq_targets.fa.gz;
KMAselect1="$REFERENCES"/"$SAMPLE".kma.select1.fasta;
KMAselect2="$REFERENCES"/"$SAMPLE".kma.select2.fasta;


# kma output parsing
cat "$MAPPING"/"$SAMPLE"."$DBshort1".tab.res | sed 's/\ /_/g' | sed 1,7d | sort -k2,2 -nr | cut -f1 -d$'\t' | head -n"$SELECT" | cut -f1 -d'.' > "$KMALIST1";
cat "$MAPPING"/"$SAMPLE"."$DBshort2".tab.res | sed 's/\ /_/g' | sed 1,7d | sort -k2,2 -nr | cut -f1 -d$'\t' | head -n"$SELECT" | sed 's/_/\ /g' > "$KMALIST2";

# construct a new fasta file containing the amount of fasta selected using "$SELECT"
filterbyname.sh -Xmx20g in="$KMAin1" out="$KMAselect1" names="$KMALIST1" include=t ow >> "$LOG1" 2>&1;
filterbyname.sh -Xmx20g in="$KMAin2" out="$KMAselect2" names="$KMALIST2" include=t ow >> "$LOG1" 2>&1;


# all input and output bbmap files
#BBMAPLIST1="$MAPPING"/"$SAMPLE".bbmap."$DBshort1".names;
#BBMAPLIST2="$MAPPING"/"$SAMPLE".bbmap."$DBshort2".names;
#BBMAPin1=/mnt/lely_DB/minikraken2/viral_20200519/library/viral/library.fna;
#BBMAPin2=/mnt/lely_DB/Old_Shared_DATABASES/VIRUSES/VirCap/VirCapSeq_targets.fa.gz;
#BBMAPselect1="$REFERENCES"/"$SAMPLE".bbmap.select1.fasta;
#BBMAPselect2="$REFERENCES"/"$SAMPLE".bbmap.select2.fasta;


# bbmap output parsing
#cat "$MAPPING"/"$SAMPLE"."$DBshort1".RPKM.tab | sed 1,6d | sed 's/\ /_/g' | sort k7,7  | cut -f1 -d$'\t' | sed 's/_/\ /g' | head -n"$SELECT" > "$BBMAPLIST1";
#cat "$MAPPING"/"$SAMPLE"."$DBshort2".RPKM.tab | sed 1,6d | sed 's/\ /_/g' | sort k7,7  | cut -f1 -d$'\t' | sed 's/_/\ /g' | head -n"$SELECT" > "$BBMAPLIST2";
# construct a new fasta file containing the amount of fasta selected using "$SELECT"
#filterbyname.sh -Xmx20g in="$BBMAPselect1" out="$BBMAPselect1" names="$BBMAPLIST1" include=t ow >> "$LOG1" 2>&1;
#filterbyname.sh -Xmx20g in="$BBMAPselect2" out="$BBMAPselect2" names="$BBMAPLIST2" include=t ow >> "$LOG1" 2>&1;

# input files
	R1="$POLISHED"/"$SAMPLE"_R1.QTR.adapter.correct.fq.gz;
	R2="$POLISHED"/"$SAMPLE"_R2.QTR.adapter.correct.fq.gz;


# map against database 1
		bbmap.sh -Xmx48g ref="$KMAselect1" build=1 in1="$R1" in2="$R2" outm="$OUTM1" "$TAGS" covstats="$COVSTATS1" rpkm="$RPKM1" untrim=t showprogress=0 ambig=best slow=t=f sam=1.3 fast=t ow >> "$LOG1" 2>&1 ;

# optional output takes very long time for processing
		samtools view -@ "$NODES" -S "$OUTM1" -b -o "$OUTMb1";
		samtools sort -@ "$NODES" "$OUTMb1" -o "$REFERENCES"/"$SAMPLE"."$DBshort1"_sorted.bam;
		samtools index "$REFERENCES"/"$SAMPLE"."$DBshort1"_sorted.bam "$REFERENCES"/"$SAMPLE"."$DBshort1"_sorted.bai;
#
# map against database 2
		bbmap.sh -Xmx48g ref="$KMAselect2" build=1 in1="$R1" in2="$R2" outm="$OUTM2" "$TAGS" covstats="$COVSTATS2" rpkm="$RPKM" untrim=t showprogress=0 ambig=best slow=t=f sam=1.3 fast=t ow >> "$LOG1" 2>&1 ;

# optional output takes very long time for processing
		samtools view -@ "$NODES" -S "$OUTM2" -b -o "$OUTMb2";
		samtools sort -@ "$NODES" "$OUTMb2" -o "$REFERENCES"/"$SAMPLE"."$DBshort2"_sorted.bam;
		samtools index "$REFERENCES"/"$SAMPLE"."$DBshort2"_sorted.bam "$REFERENCES"/"$SAMPLE"."$DBshort2"_sorted.bai;


else


# variables
TAGS='mdtag=t nhtag=t xmtag=t amtag=t nmtag=t xstag=t';
LOG1="$LOG"/"$SAMPLE".general.log;
DBshort1='viral_20200519';
DBshort2='vircap';
DB1='/mnt/lely_DB/minikraken2/viral_20200519/library/viral/library.fna';
DB2='/mnt/lely_DB/Old_Shared_DATABASES/VIRUSES/VirCap/VirCapSeq_targets.fa.gz';

# select a dedicated number of positives from the bbmap and kma output files from the different databases  

# all input and output kma files
KMALIST1="$MAPPING"/"$SAMPLE".kma."$DBshort1".names;
KMALIST2="$MAPPING"/"$SAMPLE".kma."$DBshort2".names;
KMAin1=/mnt/lely_DB/minikraken2/viral_20200519/library/viral/library.fna;
KMAin2=/mnt/lely_DB/Old_Shared_DATABASES/VIRUSES/VirCap/VirCapSeq_targets.fa.gz;
KMAselect1="$REFERENCES"/"$SAMPLE".kma.select1.fasta;
KMAselect2="$REFERENCES"/"$SAMPLE".kma.select2.fasta;


# kma output parsing
cat "$MAPPING"/"$SAMPLE"."$short".clean."$DBshort1".tab.res | sed 's/\ /_/g' | sed 1,7d | sort -k2,2 -nr | cut -f1 -d$'\t' | head -n"$SELECT" | sed 's/_/\ /g' > "$KMALIST1";
cat "$MAPPING"/"$SAMPLE"."$short".clean."$DBshort2".tab.res | sed 's/\ /_/g' | sed 1,7d | sort -k2,2 -nr | cut -f1 -d$'\t' | head -n"$SELECT" | sed 's/_/\ /g' > "$KMALIST2";

# construct a new fasta file containing the amount of fasta selected using "$SELECT"
filterbyname.sh -Xmx20g in="$KMAin1" out="$KMAselect1" names="$KMALIST1" include=t ow >> "$LOG1" 2>&1;
filterbyname.sh -Xmx20g in="$KMAin2" out="$KMAselect2" names="$KMALIST2" include=t ow >> "$LOG1" 2>&1;


# input files
	R1="$POLISHED"/"$SAMPLE"_R1.QTR."$short".clean.adapter.correct.fq.gz;
	R2="$POLISHED"/"$SAMPLE"_R2.QTR."$short".clean.adapter.correct.fq.gz;



# map against database 1
		bbmap.sh -Xmx48g ref="$KMAselect1" build=1 in1="$R1" in2="$R2" outm="$OUTM1" "$TAGS" covstats="$COVSTATS1" rpkm="$RPKM1" untrim=t showprogress=0 ambig=best slow=t=f sam=1.3 fast=t ow >> "$LOG1" 2>&1 ;

# optional output takes very long time for processing
		samtools view -@ "$NODES" -S "$OUTM1" -b -o "$OUTMb1";
		samtools sort -@ "$NODES" "$OUTMb1" -o "$REFERENCES"/"$SAMPLE"."$DBshort1"_sorted.bam;
		samtools index "$REFERENCES"/"$SAMPLE"."$DBshort1"_sorted.bam "$REFERENCES"/"$SAMPLE"."$DBshort1"_sorted.bai;
#
# map against database 2
		bbmap.sh -Xmx48g ref="$KMAselect2" build=1 in1="$R1" in2="$R2" outm="$OUTM2" "$TAGS" covstats="$COVSTATS2" rpkm="$RPKM" untrim=t showprogress=0 ambig=best slow=t=f sam=1.3 fast=t ow >> "$LOG1" 2>&1 ;

# optional output takes very long time for processing
		samtools view -@ "$NODES" -S "$OUTM2" -b -o "$OUTMb2";
		samtools sort -@ "$NODES" "$OUTMb2" -o "$REFERENCES"/"$SAMPLE"."$DBshort2"_sorted.bam;
		samtools index "$REFERENCES"/"$SAMPLE"."$DBshort2"_sorted.bam "$REFERENCES"/"$SAMPLE"."$DBshort2"_sorted.bai;


fi



	let cnt--;	
	echo -e "$cnt samples to go!";
	echo "NEXT";
	
count0=$((count0+1));

done

exit 1


##### bbmap
#
#BBMap
#Written by Brian Bushnell, from Dec. 2010 - present
#Last modified September 15, 2022
#
#Description:  Fast and accurate splice-aware read aligner.
#Please read bbmap/docs/guides/BBMapGuide.txt for more information.
#
#To index:     bbmap.sh ref=<reference fasta>
#To map:       bbmap.sh in=<reads> out=<output sam>
#To map without writing an index:
#    bbmap.sh ref=<reference fasta> in=<reads> out=<output sam> nodisk
#
#in=stdin will accept reads from standard in, and out=stdout will write to
#standard out, but file extensions are still needed to specify the format of the
#input and output files e.g. in=stdin.fa.gz will read gzipped fasta from
#standard in; out=stdout.sam.gz will write gzipped sam.
#
#Indexing Parameters (required when building the index):
#nodisk=f                Set to true to build index in memory and write nothing
#                        to disk except output.
#ref=<file>              Specify the reference sequence.  Only do this ONCE,
#                        when building the index (unless using 'nodisk').
#build=1                 If multiple references are indexed in the same directory,
#                        each needs a unique numeric ID (unless using 'nodisk').
#k=13                    Kmer length, range 8-15.  Longer is faster but uses
#                        more memory.  Shorter is more sensitive.
#                        If indexing and mapping are done in two steps, K should
#                        be specified each time.
#path=<.>                Specify the location to write the index, if you don't
#                        want it in the current working directory.
#usemodulo=f             Throw away ~80% of kmers based on remainder modulo a
#                        number (reduces RAM by 50% and sensitivity slightly).
#                        Should be enabled both when building the index AND
#                        when mapping.
#rebuild=f               Force a rebuild of the index (ref= should be set).
#
#Input Parameters:
#build=1                 Designate index to use.  Corresponds to the number
#                        specified when building the index.
#in=<file>               Primary reads input; required parameter.
#in2=<file>              For paired reads in two files.
#interleaved=auto        True forces paired/interleaved input; false forces
#                        single-ended mapping. If not specified, interleaved
#                        status will be autodetected from read names.
#fastareadlen=500        Break up FASTA reads longer than this.  Max is 500 for
#                        BBMap and 6000 for BBMapPacBio.  Only works for FASTA
#                        input (use 'maxlen' for FASTQ input).  The default for
#                        bbmap.sh is 500, and for mapPacBio.sh is 6000.
#unpigz=f                Spawn a pigz (parallel gzip) process for faster
#                        decompression than using Java.
#                        Requires pigz to be installed.
#touppercase=t           (tuc) Convert lowercase letters in reads to upper case
#                        (otherwise they will not match the reference).
#
#Sampling Parameters:
#
#reads=-1                Set to a positive number N to only process the first N
#                        reads (or pairs), then quit.  -1 means use all reads.
#samplerate=1            Set to a number from 0 to 1 to randomly select that
#                        fraction of reads for mapping. 1 uses all reads.
#skipreads=0             Set to a number N to skip the first N reads (or pairs),
#                        then map the rest.
#
#Mapping Parameters:
#fast=f                  This flag is a macro which sets other paramters to run
#                        faster, at reduced sensitivity.  Bad for RNA-seq.
#slow=f                  This flag is a macro which sets other paramters to run
#                        slower, at greater sensitivity.  'vslow' is even slower.
#maxindel=16000          Don't look for indels longer than this. Lower is faster.
#                        Set to >=100k for RNAseq with long introns like mammals.
#strictmaxindel=f        When enabled, do not allow indels longer than 'maxindel'.
#                        By default these are not sought, but may be found anyway.
#tipsearch=100           Look this far for read-end deletions with anchors
#                        shorter than K, using brute force.
#minid=0.76              Approximate minimum alignment identity to look for.
#                        Higher is faster and less sensitive.
#minhits=1               Minimum number of seed hits required for candidate sites.
#                        Higher is faster.
#local=f                 Set to true to use local, rather than global, alignments.
#                        This will soft-clip ugly ends of poor alignments.
#perfectmode=f           Allow only perfect mappings when set to true (very fast).
#semiperfectmode=f       Allow only perfect and semiperfect (perfect except for
#                        N's in the reference) mappings.
#threads=auto            (t) Set to number of threads desired.  By default, uses
#                        all cores available.
#ambiguous=best          (ambig) Set behavior on ambiguously-mapped reads (with
#                        multiple top-scoring mapping locations).
#                            best    (use the first best site)
#                            toss    (consider unmapped)
#                            random  (select one top-scoring site randomly)
#                            all     (retain all top-scoring sites)
#samestrandpairs=f       (ssp) Specify whether paired reads should map to the
#                        same strand or opposite strands.
#requirecorrectstrand=t  (rcs) Forbid pairing of reads without correct strand
#                        orientation.  Set to false for long-mate-pair libraries.
#killbadpairs=f          (kbp) If a read pair is mapped with an inappropriate
#                        insert size or orientation, the read with the lower
#                        mapping quality is marked unmapped.
#pairedonly=f            (po) Treat unpaired reads as unmapped.  Thus they will
#                        be sent to 'outu' but not 'outm'.
#rcomp=f                 Reverse complement both reads prior to mapping (for LMP
#                        outward-facing libraries).
#rcompmate=f             Reverse complement read2 prior to mapping.
#pairlen=32000           Set max allowed distance between paired reads.
#                        (insert size)=(pairlen)+(read1 length)+(read2 length)
#rescuedist=1200         Don't try to rescue paired reads if avg. insert size
#                        greater than this.  Lower is faster.
#rescuemismatches=32     Maximum mismatches allowed in a rescued read.  Lower
#                        is faster.
#averagepairdist=100     (apd) Initial average distance between paired reads.
#                        Varies dynamically; does not need to be specified.
#deterministic=f         Run in deterministic mode.  In this case it is good
#                        to set averagepairdist.  BBMap is deterministic
#                        without this flag if using single-ended reads,
#                        or run singlethreaded.
#bandwidthratio=0        (bwr) If above zero, restrict alignment band to this
#                        fraction of read length.  Faster but less accurate.
#bandwidth=0             (bw) Set the bandwidth directly.
#                        fraction of read length.  Faster but less accurate.
#usejni=f                (jni) Do alignments faster, in C code.  Requires
#                        compiling the C code; details are in /jni/README.txt.
#maxsites2=800           Don't analyze (or print) more than this many alignments
#                        per read.
#ignorefrequentkmers=t   (ifk) Discard low-information kmers that occur often.
#excludefraction=0.03    (ef) Fraction of kmers to ignore.  For example, 0.03
#                        will ignore the most common 3% of kmers.
#greedy=t                Use a greedy algorithm to discard the least-useful
#                        kmers on a per-read basis.
#kfilter=0               If positive, potential mapping sites must have at
#                        least this many consecutive exact matches.
#
#
#Quality and Trimming Parameters:
#qin=auto                Set to 33 or 64 to specify input quality value ASCII
#                        offset. 33 is Sanger, 64 is old Solexa.
#qout=auto               Set to 33 or 64 to specify output quality value ASCII
#                        offset (only if output format is fastq).
#qtrim=f                 Quality-trim ends before mapping.  Options are:
#                        'f' (false), 'l' (left), 'r' (right), and 'lr' (both).
#untrim=f                Undo trimming after mapping.  Untrimmed bases will be
#                        soft-clipped in cigar strings.
#trimq=6                 Trim regions with average quality below this
#                        (phred algorithm).
#mintrimlength=60        (mintl) Don't trim reads to be shorter than this.
#fakefastaquality=-1     (ffq) Set to a positive number 1-50 to generate fake
#                        quality strings for fasta input reads.
#ignorebadquality=f      (ibq) Keep going, rather than crashing, if a read has
#                        out-of-range quality values.
#usequality=t            Use quality scores when determining which read kmers
#                        to use as seeds.
#minaveragequality=0     (maq) Do not map reads with average quality below this.
#maqb=0                  If positive, calculate maq from this many initial bases.
#
#Output Parameters:
#out=<file>              Write all reads to this file.
#outu=<file>             Write only unmapped reads to this file.  Does not
#                        include unmapped paired reads with a mapped mate.
#outm=<file>             Write only mapped reads to this file.  Includes
#                        unmapped paired reads with a mapped mate.
#mappedonly=f            If true, treats 'out' like 'outm'.
#bamscript=<file>        (bs) Write a shell script to <file> that will turn
#                        the sam output into a sorted, indexed bam file.
#ordered=f               Set to true to output reads in same order as input.
#                        Slower and uses more memory.
#overwrite=f             (ow) Allow process to overwrite existing files.
#secondary=f             Print secondary alignments.
#sssr=0.95               (secondarysitescoreratio) Print only secondary alignments
#                        with score of at least this fraction of primary.
#ssao=f                  (secondarysiteasambiguousonly) Only print secondary
#                        alignments for ambiguously-mapped reads.
#maxsites=5              Maximum number of total alignments to print per read.
#                        Only relevant when secondary=t.
#quickmatch=f            Generate cigar strings more quickly.
#trimreaddescriptions=f  (trd) Truncate read and ref names at the first whitespace,
#                        assuming that the remainder is a comment or description.
#ziplevel=2              (zl) Compression level for zip or gzip output.
#pigz=f                  Spawn a pigz (parallel gzip) process for faster
#                        compression than Java.  Requires pigz to be installed.
#machineout=f            Set to true to output statistics in machine-friendly
#                        'key=value' format.
#printunmappedcount=f    Print the total number of unmapped reads and bases.
#                        If input is paired, the number will be of pairs
#                        for which both reads are unmapped.
#showprogress=0          If positive, print a '.' every X reads.
#showprogress2=0         If positive, print the number of seconds since the
#                        last progress update (instead of a '.').
#renamebyinsert=f        Renames reads based on their mapped insert size.
#
#Bloom-Filtering Parameters (bloomfilter.sh is the standalone version).
#bloom=f                 Use a Bloom filter to ignore reads not sharing kmers
#                        with the reference.  This uses more memory, but speeds
#                        mapping when most reads don't match the reference.
#bloomhashes=2           Number of hash functions.
#bloomminhits=3          Number of consecutive hits to be considered matched.
#bloomk=31               Bloom filter kmer length.
#bloomserial=t           Use the serialized Bloom filter for greater loading
#                        speed, if available.  If not, generate and write one.
#
#Post-Filtering Parameters:
#idfilter=0              Independant of minid; sets exact minimum identity
#                        allowed for alignments to be printed.  Range 0 to 1.
#subfilter=-1            Ban alignments with more than this many substitutions.
#insfilter=-1            Ban alignments with more than this many insertions.
#delfilter=-1            Ban alignments with more than this many deletions.
#indelfilter=-1          Ban alignments with more than this many indels.
#editfilter=-1           Ban alignments with more than this many edits.
#inslenfilter=-1         Ban alignments with an insertion longer than this.
#dellenfilter=-1         Ban alignments with a deletion longer than this.
#nfilter=-1              Ban alignments with more than this many ns.  This
#                        includes nocall, noref, and off scaffold ends.
#
#Sam flags and settings:
#noheader=f              Disable generation of header lines.
#sam=1.4                 Set to 1.4 to write Sam version 1.4 cigar strings,
#                        with = and X, or 1.3 to use M.
#saa=t                   (secondaryalignmentasterisks) Use asterisks instead of
#                        bases for sam secondary alignments.
#cigar=t                 Set to 'f' to skip generation of cigar strings (faster).
#keepnames=f             Keep original names of paired reads, rather than
#                        ensuring both reads have the same name.
#intronlen=999999999     Set to a lower number like 10 to change 'D' to 'N' in
#                        cigar strings for deletions of at least that length.
#rgid=                   Set readgroup ID.  All other readgroup fields
#                        can be set similarly, with the flag rgXX=
#                        If you set a readgroup flag to the word 'filename',
#                        e.g. rgid=filename, the input file name will be used.
#mdtag=f                 Write MD tags.
#nhtag=f                 Write NH tags.
#xmtag=f                 Write XM tags (may only work correctly with ambig=all).
#amtag=f                 Write AM tags.
#nmtag=f                 Write NM tags.
#xstag=f                 Set to 'xs=fs', 'xs=ss', or 'xs=us' to write XS tags
#                        for RNAseq using firststrand, secondstrand, or
#                        unstranded libraries.  Needed by Cufflinks.
#                        JGI mainly uses 'firststrand'.
#stoptag=f               Write a tag indicating read stop location, prefixed by YS:i:
#lengthtag=f             Write a tag indicating (query,ref) alignment lengths,
#                        prefixed by YL:Z:
#idtag=f                 Write a tag indicating percent identity, prefixed by YI:f:
#inserttag=f             Write a tag indicating insert size, prefixed by X8:Z:
#scoretag=f              Write a tag indicating BBMap's raw score, prefixed by YR:i:
#timetag=f               Write a tag indicating this read's mapping time, prefixed by X0:i:
#boundstag=f             Write a tag indicating whether either read in the pair
#                        goes off the end of the reference, prefixed by XB:Z:
#notags=f                Turn off all optional tags.
#
#Histogram and statistics output parameters:
#scafstats=<file>        Statistics on how many reads mapped to which scaffold.
#refstats=<file>         Statistics on how many reads mapped to which reference
#                        file; only for BBSplit.
#sortscafs=t             Sort scaffolds or references by read count.
#bhist=<file>            Base composition histogram by position.
#qhist=<file>            Quality histogram by position.
#aqhist=<file>           Histogram of average read quality.
#bqhist=<file>           Quality histogram designed for box plots.
#lhist=<file>            Read length histogram.
#ihist=<file>            Write histogram of insert sizes (for paired reads).
#ehist=<file>            Errors-per-read histogram.
#qahist=<file>           Quality accuracy histogram of error rates versus
#                        quality score.
#indelhist=<file>        Indel length histogram.
#mhist=<file>            Histogram of match, sub, del, and ins rates by
#                        read location.
#gchist=<file>           Read GC content histogram.
#gcbins=100              Number gchist bins.  Set to 'auto' to use read length.
#gcpairs=t               Use average GC of paired reads.
#idhist=<file>           Histogram of read count versus percent identity.
#idbins=100              Number idhist bins.  Set to 'auto' to use read length.
#statsfile=stderr        Mapping statistics are printed here.
#
#Coverage output parameters (these may reduce speed and use more RAM):
#covstats=<file>         Per-scaffold coverage info.
#rpkm=<file>             Per-scaffold RPKM/FPKM counts.
#covhist=<file>          Histogram of # occurrences of each depth level.
#basecov=<file>          Coverage per base location.
#bincov=<file>           Print binned coverage per location (one line per X bases).
#covbinsize=1000         Set the binsize for binned coverage output.
#nzo=t                   Only print scaffolds with nonzero coverage.
#twocolumn=f             Change to true to print only ID and Avg_fold instead of
#                        all 6 columns to the 'out=' file.
#32bit=f                 Set to true if you need per-base coverage over 64k.
#strandedcov=f           Track coverage for plus and minus strand independently.
#startcov=f              Only track start positions of reads.
#secondarycov=t          Include coverage of secondary alignments.
#physcov=f               Calculate physical coverage for paired reads.
#                        This includes the unsequenced bases.
#delcoverage=t           (delcov) Count bases covered by deletions as covered.
#                        True is faster than false.
#covk=0                  If positive, calculate kmer coverage statistics.
#
#Java Parameters:
#-Xmx                    This will set Java's memory usage,
#                        overriding autodetection.
#                        -Xmx20g will specify 20 gigs of RAM, and -Xmx800m
#                        will specify 800 megs.  The max is typically 85% of
#                        physical memory.  The human genome requires around 24g,
#                        or 12g with the 'usemodulo' flag.  The index uses
#                        roughly 6 bytes per reference base.
#-eoom                   This flag will cause the process to exit if an
#                        out-of-memory exception occurs.  Requires Java 8u92+.
#-da                     Disable assertions.
#
#Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter
#any problems, or post at: http://seqanswers.com/forums/showthread.php?t=41057
#
#
#####

##### samtools
#
#
#Program: samtools (Tools for alignments in the SAM format)
#Version: 1.6 (using htslib 1.6)
#
#Usage:   samtools <command> [options]
#
#Commands:
#  -- Indexing
#     dict           create a sequence dictionary file
#     faidx          index/extract FASTA
#     index          index alignment
#
#  -- Editing
#     calmd          recalculate MD/NM tags and '=' bases
#     fixmate        fix mate information
#     reheader       replace BAM header
#     rmdup          remove PCR duplicates
#     targetcut      cut fosmid regions (for fosmid pool only)
#     addreplacerg   adds or replaces RG tags
#     markdup        mark duplicates
#
#  -- File operations
#     collate        shuffle and group alignments by name
#     cat            concatenate BAMs
#     merge          merge sorted alignments
#     mpileup        multi-way pileup
#     sort           sort alignment file
#     split          splits a file by read group
#     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
#     fastq          converts a BAM to a FASTQ
#     fasta          converts a BAM to a FASTA
#
#  -- Statistics
#     bedcov         read depth per BED region
#     depth          compute the depth
#     flagstat       simple stats
#     idxstats       BAM index stats
#     phase          phase heterozygotes
#     stats          generate stats (former bamcheck)
#
#  -- Viewing
#     flags          explain BAM flags
#     tview          text alignment viewer
#     view           SAM<->BAM<->CRAM conversion
#     depad          convert padded BAM to unpadded BAM
#
#####

