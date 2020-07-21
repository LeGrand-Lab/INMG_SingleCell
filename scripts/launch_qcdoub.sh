#!/usr/bin/env bash

for i in runDoub*; do
	chmod 755 $i
	./$i
done
