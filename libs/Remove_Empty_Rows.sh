#!/bin/bash

# Simple script to remove cols that are empty
infile='${1}'
outfile='${2}'

# processes white space or tab delimted files, and will remove lines having fewer than a handful of values
awk -F'\\s+' '(NF > 3){print $0}' ${infile} > ${outfile}
