#!/bin/bash

SAMPLES="test hamzeh protein amico2"

for sample in $SAMPLES
do
	echo "testing sample $sample.nex"
	rm test/$sample.nex
	python svg_converter.py test/$sample.svg
	python test/check_trees.py test/$sample.nex test/$sample.check.nex
	
# 	if diff test/$sample.nex test/$sample.check.nex
# 		then echo "### $sample is fine"
# 		else echo "### error in $sample"
# 	fi
done

