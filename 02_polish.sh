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

			echo -e "script: 02_polish.sh" >> "$PROGVERSIONS";
			echo -e "conda environment=$ENV" >> "$PROGVERSIONS";			
			echo -e "programm\tversion" >> "$PROGVERSIONS";
			echo -e "bbmap-suite\t$PROG" >> "$PROGVERSIONS";		
			echo -e "all scripts used will be documented in this file" >> "$PROGVERSIONS";
			echo -e "\n" >> "$PROGVERSIONS";


while getopts "w:x:y:z:a:b:c:d:e:f:h:m:r:s:t:" opt; do
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
      echo "-d was triggered! $OPTARG"
      KRAKEN="`echo $(readlink -m $OPTARG)`"
      echo $KRAKEN
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

# variables
KTRIM='r';
REF="$PWD"/adapter.fa;
QTRIM='rl';
TRIMQ=20;
dbPATH=/mnt/lely_DB/AMR-Lab/contaminants/;

echo -e "short=$short";


DB="$dbPATH"/"$short";

count0=1;
countS=$(cat samples.txt | wc -l);

while [ "$count0" -le "$countS" ]; do 

	SAMPLE=$(cat samples.txt | awk 'NR=='"$count0");
	
	C="$RAW_FASTQ"/"$SAMPLE"."$short".contaminant.interleaved.fq.gz;

# input files
		R1="$RAW_FASTQ"/"$SAMPLE"_R1.fastq.gz;
		R2="$RAW_FASTQ"/"$SAMPLE"_R2.fastq.gz;
		R1c="$TMP"/"$SAMPLE"_R1.clean.fq.gz;
		R2c="$TMP"/"$SAMPLE"_R2.clean.fq.gz;

# temp output files
			CORRECTED1="$TMP"/"$SAMPLE"_R1.correct.fq.gz;
			CORRECTED2="$TMP"/"$SAMPLE"_R2.correct.fq.gz;
			ADAPTERout1="$TMP"/"$SAMPLE"_R1.adapter.fq.gz;
			ADAPTERout2="$TMP"/"$SAMPLE"_R2.adapter.fq.gz;

# output files
				OUTPUT1="$POLISHED"/"$SAMPLE"_R1.QTR.adapter.correct.fq.gz;
				OUTPUT2="$POLISHED"/"$SAMPLE"_R2.QTR.adapter.correct.fq.gz;
				OUTPUT1c="$POLISHED"/"$SAMPLE"_R1.QTR."$short".clean.adapter.correct.fq.gz;
				OUTPUT2c="$POLISHED"/"$SAMPLE"_R2.QTR."$short".clean.adapter.correct.fq.gz;
					
				HISTout="$REPORTING"/"$SAMPLE".inserthist.txt;
				LOG1="$LOG"/"$SAMPLE".general.log;

echo -e "all logs for polishing the short read data for sample = $SAMPLE" >> "$LOG1";


if [ "$CLEANUP" == 1 ];then

# READ ERROR CORRECTION
					tadpole.sh -Xmx48g in1="$R1" in2="$R2" out1="$CORRECTED1" out2="$CORRECTED2" mode=correct ow >> "$LOG1" 2>&1;
# ADAPTER TRIM
					bbduk.sh -Xmx48g in1="$CORRECTED1" in2="$CORRECTED2" out1="$ADAPTERout1" out2="$ADAPTERout2" ktrim="$KTRIM" ref="$REF" k=13 mink=6 ignorejunk=t ow >> "$LOG1" 2>&1;
# QUALITY TRIM
					bbduk.sh -Xmx48g in1="$ADAPTERout1" in2="$ADAPTERout2" out1="$OUTPUT1" out2="$OUTPUT2" qtrim="$QTRIM" trimq="$TRIMQ" minlength=25 ow >> "$LOG1" 2>&1;
# Calc insertSize from reads
					bbmerge.sh -Xmx48g in1="$OUTPUT1" in2="$OUTPUT2" ihist="$HISTout" ow >> "$LOG1" 2>&1;

else

# CONTAMINANT cleanup
					bbmap.sh -Xmx48g in1="$R1" in2="$R2" outu1="$R1c" outu2="$R2c" outm="$C"  minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 idtag printunmappedcount kfilter=25 maxsites=1 k=14 path="$DB" trimq=10 untrim ow >> "$LOG1" 2>&1;
# READ ERROR CORRECTION
					tadpole.sh -Xmx48g in1="$R1c" in2="$R2c" out1="$CORRECTED1" out2="$CORRECTED2" mode=correct ow >> "$LOG1" 2>&1;
# ADAPTER TRIM
					bbduk.sh -Xmx48g in1="$CORRECTED1" in2="$CORRECTED2" out1="$ADAPTERout1" out2="$ADAPTERout2" ktrim="$KTRIM" ref="$REF" k=13 mink=6 ignorejunk=t ow >> "$LOG1" 2>&1;
# QUALITY TRIM
					bbduk.sh -Xmx48g in1="$ADAPTERout1" in2="$ADAPTERout2" out1="$OUTPUT1c" out2="$OUTPUT2c" qtrim="$QTRIM" trimq="$TRIMQ" minlength=25 ow >> "$LOG1" 2>&1;
# Calc insertSize from reads
					bbmerge.sh -Xmx48g in1="$OUTPUT1c" in2="$OUTPUT2c" ihist="$HISTout" ow >> "$LOG1" 2>&1;

fi


# remove temp files
rm "$TMP"/*.gz;


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

##### tadpole
#
#
#Written by Brian Bushnell
#Last modified February 3, 2021
#
#Description:  Uses kmer counts to assemble contigs, extend sequences,
#or error-correct reads.  Tadpole has no upper bound for kmer length,
#but some values are not supported.  Specifically, it allows 1-31,
#multiples of 2 from 32-62, multiples of 3 from 63-93, etc.
#Please read bbmap/docs/guides/TadpoleGuide.txt for more information.
#
#Usage:
#Assembly:     tadpole.sh in=<reads> out=<contigs>
#Extension:    tadpole.sh in=<reads> out=<extended> mode=extend
#Correction:   tadpole.sh in=<reads> out=<corrected> mode=correct
#
#Recommended parameters for optimal assembly:
#tadpole.sh in=<reads> out=<contigs> shave rinse pop k=<50-70% of read length>
#
#Extension and correction may be done simultaneously.  Error correction on
#multiple files may be done like this:
#
#tadpole.sh in=libA_r1.fq,libA_merged.fq in2=libA_r2.fq,null extra=libB_r1.fq out=ecc_libA_r1.fq,ecc_libA_merged.fq out2=ecc_libA_r2.fq,null mode=correct
#
#Extending contigs with reads could be done like this:
#
#tadpole.sh in=contigs.fa out=extended.fa el=100 er=100 mode=extend extra=reads.fq k=62
#
#
#Input parameters:
#in=<file>           Primary input file for reads to use as kmer data.
#in2=<file>          Second input file for paired data.
#extra=<file>        Extra files for use as kmer data, but not for error-
#                    correction or extension.
#reads=-1            Only process this number of reads, then quit (-1 means all).
#NOTE: in, in2, and extra may also be comma-delimited lists of files.
#
#Output parameters:
#out=<file>          Write contigs (in contig mode) or corrected/extended
#                    reads (in other modes).
#out2=<file>         Second output file for paired output.
#outd=<file>         Write discarded reads, if using junk-removal flags.
#dot=<file>          Write a contigs connectivity graph (partially implemented)
#dump=<file>         Write kmers and their counts.
#fastadump=t         Write kmers and counts as fasta versus 2-column tsv.
#mincounttodump=1    Only dump kmers with at least this depth.
#showstats=t         Print assembly statistics after writing contigs.
#
#Prefiltering parameters:
#prefilter=0         If set to a positive integer, use a countmin sketch
#                    to ignore kmers with depth of that value or lower.
#prehashes=2         Number of hashes for prefilter.
#prefiltersize=0.2   (pff) Fraction of memory to use for prefilter.
#minprobprefilter=t  (mpp) Use minprob for the prefilter.
#prepasses=1         Use this many prefiltering passes; higher be more thorough
#                    if the filter is very full.  Set to 'auto' to iteratively
#                    prefilter until the remaining kmers will fit in memory.
#onepass=f           If true, prefilter will be generated in same pass as kmer
#                    counts.  Much faster but counts will be lower, by up to
#                    prefilter's depth limit.
#filtermem=0         Allows manually specifying prefilter memory in bytes, for
#                    deterministic runs.  0 will set it automatically.
#
#Hashing parameters:
#k=31                Kmer length (1 to infinity).  Memory use increases with K.
#prealloc=t          Pre-allocate memory rather than dynamically growing;
#                    faster and more memory-efficient.  A float fraction (0-1)
#                    may be specified; default is 1.
#minprob=0.5         Ignore kmers with overall probability of correctness below this.
#minprobmain=t       (mpm) Use minprob for the primary kmer counts.
#threads=X           Spawn X worker threads; default is number of logical processors.
#buildthreads=X      Spawn X contig-building threads. If not set, defaults to the same
#                    as threads.  Setting this to 1 will make contigs deterministic.
#rcomp=t             Store and count each kmer together and its reverse-complement.
#coremask=t          All kmer extensions share the same hashcode.
#fillfast=t          Speed up kmer extension lookups.
#
#Assembly parameters:
#mincountseed=3      (mcs) Minimum kmer count to seed a new contig or begin extension.
#mincountextend=2    (mce) Minimum kmer count continue extension of a read or contig.
#                    It is recommended that mce=1 for low-depth metagenomes.
#mincountretain=0    (mincr) Discard kmers with count below this.
#maxcountretain=INF  (maxcr) Discard kmers with count above this.
#branchmult1=20      (bm1) Min ratio of 1st to 2nd-greatest path depth at high depth.
#branchmult2=3       (bm2) Min ratio of 1st to 2nd-greatest path depth at low depth.
#branchlower=3       (blc) Max value of 2nd-greatest path depth to be considered low.
#minextension=2      (mine) Do not keep contigs that did not extend at least this much.
#mincontig=auto      (minc) Do not write contigs shorter than this.
#mincoverage=1       (mincov) Do not write contigs with average coverage below this.
#maxcoverage=inf     (maxcov) Do not write contigs with average coverage above this.
#trimends=0          (trim) Trim contig ends by this much.  Trimming by K/2
#                    may yield more accurate genome size estimation.
#trimcircular=t      Trim one end of contigs ending in LOOP/LOOP by K-1,
#                    to eliminate the overlapping portion.
#contigpasses=16     Build contigs with decreasing seed depth for this many iterations.
#contigpassmult=1.7  Ratio between seed depth of two iterations.
#ownership=auto      For concurrency; do not touch.
#processcontigs=f    Explore the contig connectivity graph.
#popbubbles=t        (pop) Pop bubbles; increases contiguity.  Requires
#                    additional time and memory and forces processcontigs=t.
#
#Processing modes:
#mode=contig         contig: Make contigs from kmers.
#                    extend: Extend sequences to be longer, and optionally
#                            perform error correction.
#                    correct: Error correct only.
#                    insert: Measure insert sizes.
#                    discard: Discard low-depth reads, without error correction.
#
#Extension parameters:
#extendleft=100      (el) Extend to the left by at most this many bases.
#extendright=100     (er) Extend to the right by at most this many bases.
#ibb=t               (ignorebackbranches) Do not stop at backward branches.
#extendrollback=3    Trim a random number of bases, up to this many, on reads
#                    that extend only partially.  This prevents the creation
#                    of sharp coverage discontinuities at branches.
#
#Error-correction parameters:
#ecc=f               Error correct via kmer counts.
#reassemble=t        If ecc is enabled, use the reassemble algorithm.
#pincer=f            If ecc is enabled, use the pincer algorithm.
#tail=f              If ecc is enabled, use the tail algorithm.
#eccfull=f           If ecc is enabled, use tail over the entire read.
#aggressive=f        (aecc) Use aggressive error correction settings.
#                    Overrides some other flags like errormult1 and deadzone.
#conservative=f      (cecc) Use conservative error correction settings.
#                    Overrides some other flags like errormult1 and deadzone.
#rollback=t          Undo changes to reads that have lower coverage for
#                    any kmer after correction.
#markbadbases=0      (mbb) Any base fully covered by kmers with count below
#                    this will have its quality reduced.
#markdeltaonly=t     (mdo) Only mark bad bases adjacent to good bases.
#meo=t               (markerrorreadsonly) Only mark bad bases in reads
#                    containing errors.
#markquality=0       (mq) Set quality scores for marked bases to this.
#                    A level of 0 will also convert the base to an N.
#errormult1=16       (em1) Min ratio between kmer depths to call an error.
#errormult2=2.6      (em2) Alternate ratio between low-depth kmers.
#errorlowerconst=3   (elc) Use mult2 when the lower kmer is at most this deep.
#mincountcorrect=3   (mcc) Don't correct to kmers with count under this.
#pathsimilarityfraction=0.45(psf) Max difference ratio considered similar.
#                           Controls whether a path appears to be continuous.
#pathsimilarityconstant=3   (psc) Absolute differences below this are ignored.
#errorextensionreassemble=5 (eer) Verify this many kmers before the error as
#                           having similar depth, for reassemble.
#errorextensionpincer=5     (eep) Verify this many additional bases after the
#                           error as matching current bases, for pincer.
#errorextensiontail=9       (eet) Verify additional bases before and after
#                           the error as matching current bases, for tail.
#deadzone=0          (dz) Do not try to correct bases within this distance of
#                    read ends.
#window=12           (w) Length of window to use in reassemble mode.
#windowcount=6       (wc) If more than this many errors are found within a
#                    a window, halt correction in that direction.
#qualsum=80          (qs) If the sum of the qualities of corrected bases within
#                    a window exceeds this, halt correction in that direction.
#rbi=t               (requirebidirectional) Require agreement from both
#                    directions when correcting errors in the middle part of
#                    the read using the reassemble algorithm.
#errorpath=1         (ep) For debugging purposes.
#
#Junk-removal parameters (to only remove junk, set mode=discard):
#tossjunk=f          Remove reads that cannot be used for assembly.
#                    This means they have no kmers above depth 1 (2 for paired
#                    reads) and the outermost kmers cannot be extended.
#                    Pairs are removed only if both reads fail.
#tossdepth=-1        Remove reads containing kmers at or below this depth.
#                    Pairs are removed if either read fails.
#lowdepthfraction=0  (ldf) Require at least this fraction of kmers to be
#                    low-depth to discard a read; range 0-1. 0 still
#                    requires at least 1 low-depth kmer.
#requirebothbad=f    (rbb) Only discard pairs if both reads are low-depth.
#tossuncorrectable   (tu) Discard reads containing uncorrectable errors.
#                    Requires error-correction to be enabled.
#
#Shaving parameters:
#shave=t             Remove dead ends (aka hair).
#rinse=t             Remove bubbles.
#wash=               Set shave and rinse at the same time.
#maxshavedepth=1     (msd) Shave or rinse kmers at most this deep.
#exploredist=300     (sed) Quit after exploring this far.
#discardlength=150   (sdl) Discard shavings up to this long.
#Note: Shave and rinse can produce substantially better assemblies
#for low-depth data, but they are very slow for large metagenomes.
#
#Overlap parameters (for overlapping paired-end reads only):
#merge=f             Attempt to merge overlapping reads prior to
#                    kmer-counting, and again prior to correction.  Output
#                    will still be unmerged pairs.
#ecco=f              Error correct via overlap, but do not merge reads.
#testmerge=t         Test kmer counts around the read merge junctions.  If
#                    it appears that the merge created new errors, undo it.
#
#Java Parameters:
#-Xmx                This will set Java's memory usage, overriding autodetection.
#                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
#                    The max is typically 85% of physical memory.
#-eoom               This flag will cause the process to exit if an
#                    out-of-memory exception occurs.  Requires Java 8u92+.
#-da                 Disable assertions.
#####

##### bbduk
#
#
#
#
#Written by Brian Bushnell
#Last modified July 28, 2022
#
#Description:  Compares reads to the kmers in a reference dataset, optionally
#allowing an edit distance. Splits the reads into two outputs - those that
#match the reference, and those that don't. Can also trim (remove) the matching
#parts of the reads rather than binning the reads.
#Please read bbmap/docs/guides/BBDukGuide.txt for more information.
#
#Usage:  bbduk.sh in=<input file> out=<output file> ref=<contaminant files>
#
#Input may be stdin or a fasta or fastq file, compressed or uncompressed.
#If you pipe via stdin/stdout, please include the file type; e.g. for gzipped
#fasta input, set in=stdin.fa.gz
#
#Input parameters:
#in=<file>           Main input. in=stdin.fq will pipe from stdin.
#in2=<file>          Input for 2nd read of pairs in a different file.
#ref=<file,file>     Comma-delimited list of reference files.
#                    In addition to filenames, you may also use the keywords:
#                    adapters, artifacts, phix, lambda, pjet, mtst, kapa
#literal=<seq,seq>   Comma-delimited list of literal reference sequences.
#touppercase=f       (tuc) Change all bases upper-case.
#interleaved=auto    (int) t/f overrides interleaved autodetection.
#qin=auto            Input quality offset: 33 (Sanger), 64, or auto.
#reads=-1            If positive, quit after processing X reads or pairs.
#copyundefined=f     (cu) Process non-AGCT IUPAC reference bases by making all
#                    possible unambiguous copies.  Intended for short motifs
#                    or adapter barcodes, as time/memory use is exponential.
#samplerate=1        Set lower to only process a fraction of input reads.
#samref=<file>       Optional reference fasta for processing sam files.
#
#Output parameters:
#out=<file>          (outnonmatch) Write reads here that do not contain
#                    kmers matching the database.  'out=stdout.fq' will pipe
#                    to standard out.
#out2=<file>         (outnonmatch2) Use this to write 2nd read of pairs to a
#                    different file.
#outm=<file>         (outmatch) Write reads here that fail filters.  In default
#                    kfilter mode, this means any read with a matching kmer.
#                    In any mode, it also includes reads that fail filters such
#                    as minlength, mingc, maxgc, entropy, etc.  In other words,
#                    it includes all reads that do not go to 'out'.
#outm2=<file>        (outmatch2) Use this to write 2nd read of pairs to a
#                    different file.
#outs=<file>         (outsingle) Use this to write singleton reads whose mate
#                    was trimmed shorter than minlen.
#stats=<file>        Write statistics about which contamininants were detected.
#refstats=<file>     Write statistics on a per-reference-file basis.
#rpkm=<file>         Write RPKM for each reference sequence (for RNA-seq).
#dump=<file>         Dump kmer tables to a file, in fasta format.
#duk=<file>          Write statistics in duk's format. *DEPRECATED*
#nzo=t               Only write statistics about ref sequences with nonzero hits.
#overwrite=t         (ow) Grant permission to overwrite files.
#showspeed=t         (ss) 'f' suppresses display of processing speed.
#ziplevel=2          (zl) Compression level; 1 (min) through 9 (max).
#fastawrap=70        Length of lines in fasta output.
#qout=auto           Output quality offset: 33 (Sanger), 64, or auto.
#statscolumns=3      (cols) Number of columns for stats output, 3 or 5.
#                    5 includes base counts.
#rename=f            Rename reads to indicate which sequences they matched.
#refnames=f          Use names of reference files rather than scaffold IDs.
#trd=f               Truncate read and ref names at the first whitespace.
#ordered=f           Set to true to output reads in same order as input.
#maxbasesout=-1      If positive, quit after writing approximately this many
#                    bases to out (outu/outnonmatch).
#maxbasesoutm=-1     If positive, quit after writing approximately this many
#                    bases to outm (outmatch).
#json=f              Print to screen in json format.
#
#Histogram output parameters:
#bhist=<file>        Base composition histogram by position.
#qhist=<file>        Quality histogram by position.
#qchist=<file>       Count of bases with each quality value.
#aqhist=<file>       Histogram of average read quality.
#bqhist=<file>       Quality histogram designed for box plots.
#lhist=<file>        Read length histogram.
#phist=<file>        Polymer length histogram.
#gchist=<file>       Read GC content histogram.
#enthist=<file>      Read entropy histogram.
#ihist=<file>        Insert size histogram, for paired reads in mapped sam.
#gcbins=100          Number gchist bins.  Set to 'auto' to use read length.
#maxhistlen=6000     Set an upper bound for histogram lengths; higher uses
#                    more memory.  The default is 6000 for some histograms
#                    and 80000 for others.
#
#Histograms for mapped sam/bam files only:
#histbefore=t        Calculate histograms from reads before processing.
#ehist=<file>        Errors-per-read histogram.
#qahist=<file>       Quality accuracy histogram of error rates versus quality
#                    score.
#indelhist=<file>    Indel length histogram.
#mhist=<file>        Histogram of match, sub, del, and ins rates by position.
#idhist=<file>       Histogram of read count versus percent identity.
#idbins=100          Number idhist bins.  Set to 'auto' to use read length.
#varfile=<file>      Ignore substitution errors listed in this file when
#                    calculating error rates.  Can be generated with
#                    CallVariants.
#vcf=<file>          Ignore substitution errors listed in this VCF file
#                    when calculating error rates.
#ignorevcfindels=t   Also ignore indels listed in the VCF.
#
#Processing parameters:
#k=27                Kmer length used for finding contaminants.  Contaminants
#                    shorter than k will not be found.  k must be at least 1.
#rcomp=t             Look for reverse-complements of kmers in addition to
#                    forward kmers.
#maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to
#                    increase sensitivity in the presence of errors.
#minkmerhits=1       (mkh) Reads need at least this many matching kmers
#                    to be considered as matching the reference.
#minkmerfraction=0.0 (mkf) A reads needs at least this fraction of its total
#                    kmers to hit a ref, in order to be considered a match.
#                    If this and minkmerhits are set, the greater is used.
#mincovfraction=0.0  (mcf) A reads needs at least this fraction of its total
#                    bases to be covered by ref kmers to be considered a match.
#                    If specified, mcf overrides mkh and mkf.
#hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only).
#                    Memory use is proportional to (3*K)^hdist.
#qhdist=0            Hamming distance for query kmers; impacts speed, not memory.
#editdistance=0      (edist) Maximum edit distance from ref kmers (subs
#                    and indels).  Memory use is proportional to (8*K)^edist.
#hammingdistance2=0  (hdist2) Sets hdist for short kmers, when using mink.
#qhdist2=0           Sets qhdist for short kmers, when using mink.
#editdistance2=0     (edist2) Sets edist for short kmers, when using mink.
#forbidn=f           (fn) Forbids matching of read kmers containing N.
#                    By default, these will match a reference 'A' if
#                    hdist>0 or edist>0, to increase sensitivity.
#removeifeitherbad=t (rieb) Paired reads get sent to 'outmatch' if either is
#                    match (or either is trimmed shorter than minlen).
#                    Set to false to require both.
#trimfailures=f      Instead of discarding failed reads, trim them to 1bp.
#                    This makes the statistics a bit odd.
#findbestmatch=f     (fbm) If multiple matches, associate read with sequence
#                    sharing most kmers.  Reduces speed.
#skipr1=f            Don't do kmer-based operations on read 1.
#skipr2=f            Don't do kmer-based operations on read 2.
#ecco=f              For overlapping paired reads only.  Performs error-
#                    correction with BBMerge prior to kmer operations.
#recalibrate=f       (recal) Recalibrate quality scores.  Requires calibration
#                    matrices generated by CalcTrueQuality.
#sam=<file,file>     If recalibration is desired, and matrices have not already
#                    been generated, BBDuk will create them from the sam file.
#amino=f             Run in amino acid mode.  Some features have not been
#                    tested, but kmer-matching works fine.  Maximum k is 12.
#
#Speed and Memory parameters:
#threads=auto        (t) Set number of threads to use; default is number of
#                    logical processors.
#prealloc=f          Preallocate memory in table.  Allows faster table loading
#                    and more efficient memory usage, for a large reference.
#monitor=f           Kill this process if it crashes.  monitor=600,0.01 would
#                    kill after 600 seconds under 1% usage.
#minrskip=1          (mns) Force minimal skip interval when indexing reference
#                    kmers.  1 means use all, 2 means use every other kmer, etc.
#maxrskip=1          (mxs) Restrict maximal skip interval when indexing
#                    reference kmers. Normally all are used for scaffolds<100kb,
#                    but with longer scaffolds, up to maxrskip-1 are skipped.
#rskip=              Set both minrskip and maxrskip to the same value.
#                    If not set, rskip will vary based on sequence length.
#qskip=1             Skip query kmers to increase speed.  1 means use all.
#speed=0             Ignore this fraction of kmer space (0-15 out of 16) in both
#                    reads and reference.  Increases speed and reduces memory.
#Note: Do not use more than one of 'speed', 'qskip', and 'rskip'.
#
#Trimming/Filtering/Masking parameters:
#Note - if ktrim, kmask, and ksplit are unset, the default behavior is kfilter.
#All kmer processing modes are mutually exclusive.
#Reads only get sent to 'outm' purely based on kmer matches in kfilter mode.
#
#ktrim=f             Trim reads to remove bases matching reference kmers, plus
#                    all bases to the left or right.
#                    Values:
#                       f (don't trim),
#                       r (trim to the right),
#                       l (trim to the left)
#ktrimtips=0         Set this to a positive number to perform ktrim on both
#                    ends, examining only the outermost X bases.
#kmask=              Replace bases matching ref kmers with another symbol.
#                    Allows any non-whitespace character, and processes short
#                    kmers on both ends if mink is set.  'kmask=lc' will
#                    convert masked bases to lowercase.
#maskfullycovered=f  (mfc) Only mask bases that are fully covered by kmers.
#ksplit=f            For single-ended reads only.  Reads will be split into
#                    pairs around the kmer.  If the kmer is at the end of the
#                    read, it will be trimmed instead.  Singletons will go to
#                    out, and pairs will go to outm.  Do not use ksplit with
#                    other operations such as quality-trimming or filtering.
#mink=0              Look for shorter kmers at read tips down to this length,
#                    when k-trimming or masking.  0 means disabled.  Enabling
#                    this will disable maskmiddle.
#qtrim=f             Trim read ends to remove bases with quality below trimq.
#                    Performed AFTER looking for kmers.  Values:
#                       rl (trim both ends),
#                       f (neither end),
#                       r (right end only),
#                       l (left end only),
#                       w (sliding window).
#trimq=6             Regions with average quality BELOW this will be trimmed,
#                    if qtrim is set to something other than f.  Can be a
#                    floating-point number like 7.3.
#trimclip=f          Trim soft-clipped bases from sam files.
#minlength=10        (ml) Reads shorter than this after trimming will be
#                    discarded.  Pairs will be discarded if both are shorter.
#mlf=0               (minlengthfraction) Reads shorter than this fraction of
#                    original length after trimming will be discarded.
#maxlength=          Reads longer than this after trimming will be discarded.
#minavgquality=0     (maq) Reads with average quality (after trimming) below
#                    this will be discarded.
#maqb=0              If positive, calculate maq from this many initial bases.
#minbasequality=0    (mbq) Reads with any base below this quality (after
#                    trimming) will be discarded.
#maxns=-1            If non-negative, reads with more Ns than this
#                    (after trimming) will be discarded.
#mcb=0               (minconsecutivebases) Discard reads without at least
#                    this many consecutive called bases.
#ottm=f              (outputtrimmedtomatch) Output reads trimmed to shorter
#                    than minlength to outm rather than discarding.
#tp=0                (trimpad) Trim this much extra around matching kmers.
#tbo=f               (trimbyoverlap) Trim adapters based on where paired
#                    reads overlap.
#strictoverlap=t     Adjust sensitivity for trimbyoverlap mode.
#minoverlap=14       Require this many bases of overlap for detection.
#mininsert=40        Require insert size of at least this for overlap.
#                    Should be reduced to 16 for small RNA sequencing.
#tpe=f               (trimpairsevenly) When kmer right-trimming, trim both
#                    reads to the minimum length of either.
#forcetrimleft=0     (ftl) If positive, trim bases to the left of this position
#                    (exclusive, 0-based).
#forcetrimright=0    (ftr) If positive, trim bases to the right of this position
#                    (exclusive, 0-based).
#forcetrimright2=0   (ftr2) If positive, trim this many bases on the right end.
#forcetrimmod=0      (ftm) If positive, right-trim length to be equal to zero,
#                    modulo this number.
#restrictleft=0      If positive, only look for kmer matches in the
#                    leftmost X bases.
#restrictright=0     If positive, only look for kmer matches in the
#                    rightmost X bases.
#NOTE:  restrictleft and restrictright are mutually exclusive.  If trimming
#       both ends is desired, use ktrimtips.
#mingc=0             Discard reads with GC content below this.
#maxgc=1             Discard reads with GC content above this.
#gcpairs=t           Use average GC of paired reads.
#                    Also affects gchist.
#tossjunk=f          Discard reads with invalid characters as bases.
#swift=f             Trim Swift sequences: Trailing C/T/N R1, leading G/A/N R2.
#
#Header-parsing parameters - these require Illumina headers:
#chastityfilter=f    (cf) Discard reads with id containing ' 1:Y:' or ' 2:Y:'.
#barcodefilter=f     Remove reads with unexpected barcodes if barcodes is set,
#                    or barcodes containing 'N' otherwise.  A barcode must be
#                    the last part of the read header.  Values:
#                       t:     Remove reads with bad barcodes.
#                       f:     Ignore barcodes.
#                       crash: Crash upon encountering bad barcodes.
#barcodes=           Comma-delimited list of barcodes or files of barcodes.
#xmin=-1             If positive, discard reads with a lesser X coordinate.
#ymin=-1             If positive, discard reads with a lesser Y coordinate.
#xmax=-1             If positive, discard reads with a greater X coordinate.
#ymax=-1             If positive, discard reads with a greater Y coordinate.
#
#Polymer trimming:
#trimpolya=0         If greater than 0, trim poly-A or poly-T tails of
#                    at least this length on either end of reads.
#trimpolygleft=0     If greater than 0, trim poly-G prefixes of at least this
#                    length on the left end of reads.  Does not trim poly-C.
#trimpolygright=0    If greater than 0, trim poly-G tails of at least this
#                    length on the right end of reads.  Does not trim poly-C.
#trimpolyg=0         This sets both left and right at once.
#filterpolyg=0       If greater than 0, remove reads with a poly-G prefix of
#                    at least this length (on the left).
#Note: there are also equivalent poly-C flags.
#
#Polymer tracking:
#pratio=base,base    'pratio=G,C' will print the ratio of G to C polymers.
#plen=20             Length of homopolymers to count.
#
#Entropy/Complexity parameters:
#entropy=-1          Set between 0 and 1 to filter reads with entropy below
#                    that value.  Higher is more stringent.
#entropywindow=50    Calculate entropy using a sliding window of this length.
#entropyk=5          Calculate entropy using kmers of this length.
#minbasefrequency=0  Discard reads with a minimum base frequency below this.
#entropytrim=f       Values:
#                       f:  (false) Do not entropy-trim.
#                       r:  (right) Trim low entropy on the right end only.
#                       l:  (left) Trim low entropy on the left end only.
#                       rl: (both) Trim low entropy on both ends.
#entropymask=f       Values:
#                       f:  (filter) Discard low-entropy sequences.
#                       t:  (true) Mask low-entropy parts of sequences with N.
#                       lc: Change low-entropy parts of sequences to lowercase.
#entropymark=f       Mark each base with its entropy value.  This is on a scale
#                    of 0-41 and is reported as quality scores, so the output
#                    should be fastq or fasta+qual.
#NOTE: If set, entropytrim overrides entropymask.
#
#Cardinality estimation:
#cardinality=f       (loglog) Count unique kmers using the LogLog algorithm.
#cardinalityout=f    (loglogout) Count unique kmers in output reads.
#loglogk=31          Use this kmer length for counting.
#loglogbuckets=2048  Use this many buckets for counting.
#khist=<file>        Kmer frequency histogram; plots number of kmers versus
#                    kmer depth.  This is approximate.
#khistout=<file>     Kmer frequency histogram for output reads.
#
#Java Parameters:
#
#-Xmx                This will set Java's memory usage, overriding autodetection.
#                    -Xmx20g will
#                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
#                    The max is typically 85% of physical memory.
#-eoom               This flag will cause the process to exit if an
#                    out-of-memory exception occurs.  Requires Java 8u92+.
#-da                 Disable assertions.
#
#Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
#
#
#####

##### bbmerge
#
#
#
#Written by Brian Bushnell and Jonathan Rood
#Last modified June 26, 2019
#
#Description:  Merges paired reads into single reads by overlap detection.
#With sufficient coverage, can merge nonoverlapping reads by kmer extension.
#Kmer modes (Tadpole or Bloom Filter) require much more memory, and should
#be used with the bbmerge-auto.sh script rather than bbmerge.sh.
#Please read bbmap/docs/guides/BBMergeGuide.txt for more information.
#
#Usage for interleaved files:    bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>
#Usage for paired files:         bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>
#
#Input may be stdin or a file, fasta or fastq, raw or gzipped.
#
#Input parameters:
#in=null              Primary input. 'in2' will specify a second file.
#interleaved=auto     May be set to true or false to override autodetection of
#                     whether the input file as interleaved.
#reads=-1             Quit after this many read pairs (-1 means all).
#
#Output parameters:
#out=<file>           File for merged reads. 'out2' will specify a second file.
#outu=<file>          File for unmerged reads. 'outu2' will specify a second file.
#outinsert=<file>     (outi) File to write read names and insert sizes.
#outadapter=<file>    (outa) File to write consensus adapter sequences.
#outc=<file>          File to write input read kmer cardinality estimate.
#ihist=<file>         (hist) Insert length histogram output file.
#nzo=t                Only print histogram bins with nonzero values.
#showhiststats=t      Print extra header lines with statistical information.
#ziplevel=2           Set to 1 (lowest) through 9 (max) to change compression
#                     level; lower compression is faster.
#ordered=f            Output reads in same order as input.
#mix=f                Output both the merged (or mergable) and unmerged reads
#                     in the same file (out=).  Useful for ecco mode.
#
#Trimming/Filtering parameters:
#qtrim=f              Trim read ends to remove bases with quality below minq.
#                     Trims BEFORE merging.
#                     Values: t (trim both ends),
#                             f (neither end),
#                             r (right end only),
#                             l (left end only).
#qtrim2=f             May be specified instead of qtrim to perform trimming
#                     only if merging is unsuccessful, then retry merging.
#trimq=10             Trim quality threshold.  This may be a comma-delimited
#                     list (ascending) to try multiple values.
#minlength=1          (ml) Reads shorter than this after trimming, but before
#                     merging, will be discarded. Pairs will be discarded only
#                     if both are shorter.
#maxlength=-1         Reads with longer insert sizes will be discarded.
#tbo=f                (trimbyoverlap) Trim overlapping reads to remove
#                     rightmost (3') non-overlapping portion, instead of joining.
#minavgquality=0      (maq) Reads with average quality below this, after
#                     trimming, will not be attempted to be merged.
#maxexpectederrors=0  (mee) If positive, reads with more combined expected
#                     errors than this will not be attempted to be merged.
#forcetrimleft=0      (ftl) If nonzero, trim left bases of the read to
#                     this position (exclusive, 0-based).
#forcetrimright=0     (ftr) If nonzero, trim right bases of the read
#                     after this position (exclusive, 0-based).
#forcetrimright2=0    (ftr2) If positive, trim this many bases on the right end.
#forcetrimmod=5       (ftm) If positive, trim length to be equal to
#                     zero modulo this number.
#ooi=f                Output only incorrectly merged reads, for testing.
#trimpolya=t          Trim trailing poly-A tail from adapter output.  Only
#                     affects outadapter.  This also trims poly-A followed
#                     by poly-G, which occurs on NextSeq.
#
#Processing Parameters:
#usejni=f             (jni) Do overlapping in C code, which is faster.  Requires
#                     compiling the C code; details are in /jni/README.txt.
#                     However, the jni path is currently disabled.
#merge=t              Create merged reads.  If set to false, you can still
#                     generate an insert histogram.
#ecco=f               Error-correct the overlapping part, but don't merge.
#trimnonoverlapping=f (tno) Trim all non-overlapping portions, leaving only
#                     consensus sequence.  By default, only sequence to the
#                     right of the overlap (adapter sequence) is trimmed.
#useoverlap=t         Attempt find the insert size using read overlap.
#mininsert=35         Minimum insert size to merge reads.
#mininsert0=35        Insert sizes less than this will not be considered.
#                     Must be less than or equal to mininsert.
#minoverlap=12        Minimum number of overlapping bases to allow merging.
#minoverlap0=8        Overlaps shorter than this will not be considered.
#                     Must be less than or equal to minoverlap.
#minq=9               Ignore bases with quality below this.
#maxq=41              Cap output quality scores at this.
#entropy=t            Increase the minimum overlap requirement for low-
#                     complexity reads.
#efilter=6            Ban overlaps with over this many times the expected
#                     number of errors.  Lower is more strict. -1 disables.
#pfilter=0.00004      Ban improbable overlaps.  Higher is more strict. 0 will
#                     disable the filter; 1 will allow only perfect overlaps.
#kfilter=0            Ban overlaps that create kmers with count below
#                     this value (0 disables).  If this is used minprob should
#                     probably be set to 0.  Requires good coverage.
#ouq=f                Calculate best overlap using quality values.
#owq=t                Calculate best overlap without using quality values.
#usequality=t         If disabled, quality values are completely ignored,
#                     both for overlap detection and filtering.  May be useful
#                     for data with inaccurate quality values.
#iupacton=f           (itn) Change ambiguous IUPAC symbols to N.
#adapter=             Specify the adapter sequences used for these reads, if
#                     known; this can be a fasta file or a literal sequence.
#                     Read 1 and 2 can have adapters specified independently
#                     with the adapter1 and adapter2 flags.  adapter=default
#                     will use a list of common adapter sequences.
#
#Ratio Mode:
#ratiomode=t          Score overlaps based on the ratio of matching to
#                     mismatching bases.
#maxratio=0.09        Max error rate; higher increases merge rate.
#ratiomargin=5.5      Lower increases merge rate; min is 1.
#ratiooffset=0.55     Lower increases merge rate; min is 0.
#maxmismatches=20     Maximum mismatches allowed in overlapping region.
#ratiominoverlapreduction=3  This is the difference between minoverlap in
#                     flat mode and minoverlap in ratio mode; generally,
#                     minoverlap should be lower in ratio mode.
#minsecondratio=0.1   Cutoff for second-best overlap ratio.
#forcemerge=f         Disable all filters and just merge everything
#                     (not recommended).
#
#Flat Mode:
#flatmode=f           Score overlaps based on the total number of mismatching
#                     bases only.
#margin=2             The best overlap must have at least 'margin' fewer
#                     mismatches than the second best.
#mismatches=3         Do not allow more than this many mismatches.
#requireratiomatch=f  (rrm) Require the answer from flat mode and ratio mode
#                     to agree, reducing false positives if both are enabled.
#trimonfailure=t      (tof) If detecting insert size by overlap fails,
#                     the reads will be trimmed and this will be re-attempted.
#
#
#*** Ratio Mode and Flat Mode may be used alone or simultaneously. ***
#*** Ratio Mode is usually more accurate and is the default mode. ***
#
#
#Strictness (these are mutually exclusive macros that set other parameters):
#strict=f             Decrease false positive rate and merging rate.
#verystrict=f         (vstrict) Greatly decrease FP and merging rate.
#ultrastrict=f        (ustrict) Decrease FP and merging rate even more.
#maxstrict=f          (xstrict) Maximally decrease FP and merging rate.
#loose=f              Increase false positive rate and merging rate.
#veryloose=f          (vloose) Greatly increase FP and merging rate.
#ultraloose=f         (uloose) Increase FP and merging rate even more.
#maxloose=f           (xloose) Maximally decrease FP and merging rate.
#fast=f               Fastest possible mode; less accurate.
#
#Tadpole Parameters (for read extension and error-correction):
#*Note: These require more memory and should be run with bbmerge-auto.sh.*
#k=31                 Kmer length.  31 (or less) is fastest and uses the least
#                     memory, but higher values may be more accurate.
#                     60 tends to work well for 150bp reads.
#extend=0             Extend reads to the right this much before merging.
#                     Requires sufficient (>5x) kmer coverage.
#extend2=0            Extend reads this much only after a failed merge attempt,
#                     or in rem/rsem mode.
#iterations=1         (ei) Iteratively attempt to extend by extend2 distance
#                     and merge up to this many times.
#rem=f                (requireextensionmatch) Do not merge if the predicted
#                     insert size differs before and after extension.
#                     However, if only the extended reads overlap, then that
#                     insert will be used.  Requires setting extend2.
#rsem=f               (requirestrictextensionmatch) Similar to rem but stricter.
#                     Reads will only merge if the predicted insert size before
#                     and after extension match.  Requires setting extend2.
#                     Enables the lowest possible false-positive rate.
#ecctadpole=f         (ecct) If reads fail to merge, error-correct with Tadpole
#                     and try again.  This happens prior to extend2.
#reassemble=t         If ecct is enabled, use Tadpole's reassemble mode for
#                     error correction.  Alternatives are pincer and tail.
#removedeadends       (shave) Remove kmers leading to dead ends.
#removebubbles        (rinse) Remove kmers in error bubbles.
#mindepthseed=3       (mds) Minimum kmer depth to begin extension.
#mindepthextend=2     (mde) Minimum kmer depth continue extension.
#branchmult1=20       Min ratio of 1st to 2nd-greatest path depth at high depth.
#branchmult2=3        Min ratio of 1st to 2nd-greatest path depth at low depth.
#branchlower=3        Max value of 2nd-greatest path depth to be considered low.
#ibb=t                Ignore backward branches when extending.
#extra=<file>         A file or comma-delimited list of files of reads to use
#                     for kmer counting, but not for merging or output.
#prealloc=f           Pre-allocate memory rather than dynamically growing;
#                     faster and more memory-efficient for large datasets.
#                     A float fraction (0-1) may be specified, default 1.
#prefilter=0          If set to a positive integer, use a countmin sketch to
#                     ignore kmers with depth of that value or lower, to
#                     reduce memory usage.
#filtermem=0          Allows manually specifying prefilter memory in bytes, for
#                     deterministic runs.  0 will set it automatically.
#minprob=0.5          Ignore kmers with overall probability of correctness
#                     below this, to reduce memory usage.
#minapproxoverlap=26  For rem mode, do not merge reads if the extended reads
#                     indicate that the raw reads should have overlapped by
#                     at least this much, but no overlap was found.
#
#
#Bloom Filter Parameters (for kmer operations with less memory than Tadpole)
#*Note: These require more memory and should be run with bbmerge-auto.sh.*
#eccbloom=f           (eccb) If reads fail to merge, error-correct with bbcms
#                     and try again.
#testmerge=f          Test kmer counts around the read merge junctions.  If
#                     it appears that the merge created new errors, undo it.
#                     This reduces the false-positive rate, but not as much as
#                     rem or rsem.
#
#Java Parameters:
#-Xmx                 This will set Java's memory usage,
#                     overriding autodetection.
#                     For example, -Xmx400m will specify 400 MB RAM.
#-eoom                This flag will cause the process to exit if an
#                     out-of-memory exception occurs.  Requires Java 8u92+.
#-da                  Disable assertions.
#
#Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
#
#
#####



