#!/usr/bin/env bash

printf "**§**\n This is a MODULAR QC, it uses the scripts:\n\
	runQCetDoub_splittedSCRAN.R\n\
	runDoub_splittedFINDER.R\n\
    BEFORE RUNNING: set files paths and input variables directly into these two R scripts.\n\
    CAUTION: save previous work otherwise it will be overwritted!! \n **§** \n\n"\


for i in run*Doub*; do
	chmod 755 $i
	echo "Running" $i
	#./$i
done
