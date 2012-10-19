#!/bin/bash
# version 2 
#
# runcode, runs program x times inserting new parameters
# and renaming and moving output data files

. basicheader
declare PROJECT="diffusion"
declare OUTFILE="phys06a.dat"
declare INITIALS="0.01"                #initial conditions that do not change over loop
PARAM=(1 2 3)

# This loop does the moving and renaming
for((COUNTER=0; COUNTER<3; COUNTER++)); do
./$PROJECT $INITIALS ${PARAM[COUNTER]}
echo "Running $PROJECT for $INITIALS and ${PARAM[COUNTER]}"
mv $OUTFILE $DATADIR/"$PROJECT"_"$INITIALS"_"${PARAM[COUNTER]}.dat"
done



exit 0