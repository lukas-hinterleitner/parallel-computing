#!/bin/bash

if ! [[ "$1" == "omp" || "$1" == "mpi" ]]; then
  echo "argument is missing. it is either 'atomic', 'mutex' or 'custom'";
  exit 1;
fi

echo "start compiling...";

if [[ $1 == "omp" ]]; then
  /opt/global/gcc-11.2.0/bin/g++ a2-omp.cpp -O2 -std=c++20 -fopenmp -lm -o a2-omp || exit 1;
#else
  # /opt/global/gcc-11.2.0/bin/g++ -O2 -lpthread -std=c++20 -o a2-"$1" a2-"$1".cpp || exit 1;
  # todo: add MPI compilation
fi

echo "compiling finished";
echo "=============================================="

echo "measuring execution time for 1 - 32 threads with '$1'";

truncate -s 0 execution-time-"$1".csv;

counter=1;

echo "num_thread,execution-time" >> execution-time-"$1".csv;

while [ $counter -le 32 ]
do
  export OMP_NUM_THREADS=$counter;
	time=$(srun --nodes=1 ./a2-"$1" --num-threads $counter --no-verify);
	echo "$time";
	echo "${time//$'\t'/','}" >> execution-time-"$1".csv;

	((counter++));
done
