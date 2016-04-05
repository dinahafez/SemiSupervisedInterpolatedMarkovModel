#!/usr/bin/bash



##Finding out the path
sPath="`dirname \"$0\"`"
sPath="`( cd \"$sPath\" && pwd )`"

wdir=$(mktemp -d)


peakfile=$1
regionfile=$2
outdir=$3

echo "sorting files..."
bedtools sort -i $peakfile > $wdir/peakfile.bed
bedtools sort -i $regionfile > $wdir/regionfile.bed

echo "Getting intersect..."
bedtools intersect -wb -a $wdir/peakfile.bed -b $wdir/regionfile.bed > $wdir/intersect.bed

echo "Getting results..."
mkdir $outdir
Rscript $sPath/peakcov.r -intersect=$wdir/intersect.bed -out=$outdir
