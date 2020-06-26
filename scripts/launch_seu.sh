#!/usr/bin/bash
pwd
for i in $(ls seu*);do
	chmod 755 ${i}
done

for i in $(ls seu*);do
	echo """=== ${i} ====="
	./${i}

done
