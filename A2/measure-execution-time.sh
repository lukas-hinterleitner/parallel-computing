#!/bin/bash

# check if argument 1 is either 'omp' or 'mpi' and check if argument 2 and 3 are numbers
if ! [[ "$1" == "omp" ]] && [[ "$2" =~ ^[0-9]+$ && "$3" =~ ^[0-9]+$ && "$4" =~ ^[0-9]+$ ]]; then
  echo "argument is missing. it is either 'omp n m iteration'";
  exit 1;
fi

echo "start compiling...";

if [[ $1 == "omp" ]]; then
  /opt/global/gcc-11.2.0/bin/g++ a2-omp.cpp -O2 -std=c++20 -fopenmp -lm -o a2-omp || exit 1;
fi

epsilon=0.01;

echo "compiling finished";
echo "=============================================="

echo "measuring execution time for 1 - 32 threads with '$1' --n '$2' --m '$3' --epsilon '$epsilon' --max-iterations '$4''";

filename=execution-time-"$1"-"$2"-"$3"-"$4"-"$epsilon".csv;
truncate -s 0 "$filename";

counter=1;

echo "num_thread,execution-time,iterations" >> "$filename";

export OMP_NESTED=true;

while [ $counter -le 32 ]
do
  export OMP_NUM_THREADS=$counter;
  echo "num threads: $counter";
  time=$(srun --nodes=1 ./a2-"$1" --n "$2" --m "$3" --epsilon "$epsilon" --max-iterations "$4" --no-verify);
  echo "$time";
  IFS=' ';
  read -ra strarr <<< "$time";
  echo "$counter,${strarr[2]},${strarr[5]}" >> "$filename";
  #echo "${time//$'\t'/','}" >> execution-time-"$1".csv;

  ((counter++));
done


