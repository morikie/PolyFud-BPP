#!/bin/bash

#If you are behind a proxy add the following lines to your ~/.wgetrc:
#use_proxy=on
#ftp_proxy=proxy_address:port


declare -a chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

ucscFtp="ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
zipExtension=".fa.gz"
extension=".fa"
refGen="reference_genome.fa"
targetFolder="../bin/reference_genome/hg19/"

for i in "${chromosomes[@]}" 
do
	url=$ucscFtp$i$zipExtension
	zippedFile=$i$zipExtension
	if [ ! -f $i$extension ]; then
		wget -v $url -P $targetFolder
		gunzip $targetFolder$zippedFile
	fi
	tr '[:upper:]' '[:lower:]' < $targetFolder$i$extension > $targetFolder$i$extension'temp'
	mv $targetFolder$i$extension'temp' $targetFolder$i$extension
	cat "$targetFolder$i$extension" >> $targetFolder$refGen 
done


