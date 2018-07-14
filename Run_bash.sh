#!/bin/bash


let CMIN=22
let CMAX=25
#let CMIN=0
#let CMAX=93

# We will run the python thing from cloud CMIN-1 to cloud CMAX-1
COUNTER=$CMIN
while [ $COUNTER -le $CMAX ]; do
	ipython Run_bash.py $COUNTER
	let COUNTER=COUNTER+1
done

#ipython Run_bash.py 91
#ipython Run_bash.py 104
