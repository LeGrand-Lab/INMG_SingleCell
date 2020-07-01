#!/usr/bin/env bash
# Launch seurat+FItSNE on homeostasis data (D0)
# each dataset independently
# assumes compressed folders in 'data/'
# Runs from 'scripts/'
# Joha GL 2020

pwd
MYPATH=$(pwd)  # expected home/MYHOME/scripts

echo "uncompress data keeping .tar.gz files, move to data"
cd ../data
for i in $(ls *D0.tar.gz); do
	tar -zxvf $i 
	echo $i
done	

echo "return to scripts/"
cd $MYPATH
for i in $(ls seu*D0.R);do
	chmod 755 ${i}
	echo $i
done

echo "execute D0 seurat analysis +FItSNE visuals" 
for i in $(ls seu*D0.R);do
	echo """=== ${i} ====="
	./${i}
done

echo "cd data again"
cd ../data
for d in */ ;do
	if [ ! -f ${d%/}".tar.gz" ]; then
		tar -czvf ${d%/}".tar.gz" $d
	fi
done

echo "please remember NOT TO commit any voluminous uncompressed content from 'data/'"

echo "END"	


