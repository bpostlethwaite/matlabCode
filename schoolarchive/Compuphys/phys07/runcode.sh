#!/bin/bash
# version 2 
#
# runcode, runs program x times inserting new parameters
# and renaming and moving output data files

. basicheader
declare PROJECT="propane"
declare OUTFILE="phys07.dat"
declare INITIALS="1"                #initial conditions that do not change over loop
PARAM=(1)

# This loop does the moving and renaming
for((COUNTER=0; COUNTER<1; COUNTER++)); do
./$PROJECT $INITIALS ${PARAM[COUNTER]}
echo "Running $PROJECT for $INITIALS and ${PARAM[COUNTER]}"
mv $OUTFILE $DATADIR/"$PROJECT"_"$INITIALS"_"${PARAM[COUNTER]}.dat"
#mv $OUTFILE $DATADIR/ConstantHeatTest1.dat
done



exit 0