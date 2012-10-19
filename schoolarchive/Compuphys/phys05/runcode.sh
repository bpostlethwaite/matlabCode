#!/bin/bash
# version 2 
#
# runcode, runs program x times inserting new parameters
# and renaming and moving output data files

declare PROJECT="driven"
declare OUTFILE="phys05.dat"
declare INITIALS="0.2 0"                #initial conditions that do not change over loop
PARAM=(1.35 1.44 1.465)

# This loop does the moving and renaming
for((COUNTER=0; COUNTER<3; COUNTER++)); do
./$PROJECT $INITIALS ${PARAM[COUNTER]}
echo "${PARAM[COUNTER]}"
mv $OUTFILE $DATADIR/$PROJECT${PARAM[COUNTER]}.dat
done

#gnuplot -persist driven1.gp

exit 0