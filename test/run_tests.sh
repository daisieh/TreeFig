#!/bin/bash

SAMPLES="test hamzeh protein amico2"

for sample in $SAMPLES
do
	rm test/$sample.nex
	python svg_converter.py test/$sample.svg
	if diff test/$sample.nex test/$sample.check.nex
		then echo "### $sample is fine"
		else echo "### error in $sample"
	fi
done

