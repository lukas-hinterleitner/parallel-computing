#!/bin/bash

if ! [[ "$1" == "atomic" || "$1" == "mutex" || "$1" == "custom" ]]; then
  echo "argument is missing. it is either 'atomic', 'mutex' or 'custom'";
  exit 1;
fi

echo "start compiling..";

/opt/global/gcc-11.2.0/bin/g++ -O2 -lpthread -std=c++20 -o a1-parallel-"$1" a1-parallel-"$1".cpp || echo "compiling failed"; exit 1;

echo "compiling finished";
echo "=============================================="

echo "measuring execution time for 1 - 32 threads with '$1'";

truncate -s 0 execution-time-"$1".csv;

counter=1;

echo "num_thread,execution-time" >> execution-time-"$1".csv;

while [ $counter -le 32 ]
do
	time=$(srun --nodes=1 ./a1-parallel-"$1" --num-threads $counter --only-exec-times);
	echo "$time";
	echo "${time//$'\t'/','}" >> execution-time-"$1".csv;

	((counter++));
done
