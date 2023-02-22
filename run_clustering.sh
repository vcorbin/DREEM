#!/usr/bin/env bash
FASTQ=$1
REF=$2
PRIMERS=$3

# Bowtie2 command
BAMFILE="${FASTQ%.*}".bam
if [ ! -f $BAMFILE ] ; then
    # alignment command to temp_bamfile echo $FASTQ > temp_$BAMFILE
    mv temp_$BAMFILE $BAMFILE
else
    echo "Skipping the alignment"
fi

# Make bitvector
python make_bitvector.py $BAMFILE -r $2 -p $3

# Run clustering
parallel.py