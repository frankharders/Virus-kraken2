#!/bin/bash

while getopts "w:x:a:b:c:d:e:f:m:r:s:t:" opt; do
  case $opt in
     w)
      echo "-w was triggered! $OPTARG"
      WORKDIR="`echo $(readlink -m $OPTARG)`"
      echo $WORKDIR
      ;;
	 x)
      echo "-x was triggered! $OPTARG"
      CLEANUP="`echo $(readlink -m $OPTARG)`"
      echo $CLEANUP
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
      echo "-f was triggered! $OPTARG"
      REFERENCES="`echo $(readlink -m $OPTARG)`"
      echo $REFERENCES
      ;;	  
	 m)
      echo "-m was triggered! $OPTARG"
      MAPPING="`echo $(readlink -m $OPTARG)`"
      echo $MAPPING
      ;;
	 r)
      echo "-r was triggered! $OPTARG"
      REPORTING="`echo $(readlink -m $OPTARG)`"
      echo $REPORTING
      ;;
     s)
      echo "-r was triggered! $OPTARG"
      TMP="`echo $(readlink -m $OPTARG)`"
      echo $TMP
      ;;
     t)
      echo "-r was triggered! $OPTARG"
      LOG="`echo $(readlink -m $OPTARG)`"
      echo $LOG
      ;;	  
	\?)
      echo "$OPTARG" >&2
      ;;

  esac
done

#create several directories

mkdir -p "$RAW_FASTQ" "$RAWSTATS" "$POLISHED" "$TRIMMEDSTATS" "$KRAKEN" "$REFERENCES" "$MAPPING" "$REPORTING" "$TMP" "$LOG";


exit 1

