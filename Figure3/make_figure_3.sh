#!/bin/bash

declare -a codes_=('Code/LV.py' 'Code/FC.py' 'Code/CR.py' 'Code/HS.py' 'Code/ML.py' 'Code/A.py' 'Code/K.py' 'Code/D.py')

tmp=0
### number of realization of the process
number_of_iterations=1000
echo 'Running the statistics for each model' $number_of_iterations 'times'
for j in "${codes_[@]}"
do

	python -W ignore $j $number_of_iterations -o 3 > /dev/null 2>&1
	tmp=$((tmp + 1))
	echo 'model' $tmp':' $j ' --- done'
done

echo 'Now making the Figure...'
mv *txt Figure/
cd Figure/
chmod +x mk_fg.sh
./mk_fg.sh
rm *txt *aux *log
