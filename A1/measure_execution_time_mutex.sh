#!/bin/bash

truncate -s 0 execution_time_mutex.csv

counter=1;

echo "num_thread,execution_time" >> execution_time_mutex.csv

while [ $counter -le 32 ]
do
	time=$(srun --nodes=1 ./a1-parallel-atomic --num-threads $counter --only-exec-times)
	echo "$time";
	echo "${time//$'\t'/','}" >> execution_time_mutex.csv

	((counter++))
done
