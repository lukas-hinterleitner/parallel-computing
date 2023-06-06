#!/bin/bash

# check if argument 1 is either 'omp' or 'mpi' and check if argument 2 and 3 are numbers
if ! [[ "$1" =~ ^[0-9]+$ && "$2" =~ ^[0-9]+$ && "$3" =~ ^[0-9]+$ ]]; then
  echo "arguments are missing. they are 'm n iteration'";
  exit 1;
fi

echo "start compiling...";

/opt/global/gcc-11.2.0/bin/g++ a2-omp.cpp -O2 -std=c++20 -fopenmp -lm -o a2-omp || exit 1;

epsilon=0.01;

echo "compiling finished";
echo "=============================================="

echo "measuring execution time for 1 - 32 threads with omp --m '$1' --n '$2' --epsilon '$epsilon' --max-iterations '$3''";

filename=execution-time-omp-dynamic-"$1"-"$2"-"$3"-"$epsilon".csv;
truncate -s 0 "$filename";

counter=1;

echo "num_thread,execution-time,iterations" >> "$filename";

export OMP_NESTED=true;

export OMP_SCHEDULE=dynamic;

while [ $counter -le 32 ]
do
  export OMP_NUM_THREADS=$counter;
  echo "num threads: $counter";
  time=$(srun --nodes=1 ./a2-omp --m "$1" --n "$2" --epsilon "$epsilon" --max-iterations "$3" --no-verify);
  echo "$time";
  IFS=' ';
  read -ra strarr <<< "$time";
  echo "$counter,${strarr[2]},${strarr[5]}" >> "$filename";
  #echo "${time//$'\t'/','}" >> execution-time-"$1".csv;

  ((counter++));
done

echo "starting static computation";

filename=execution-time-omp-static-"$1"-"$2"-"$3"-"$epsilon".csv;
truncate -s 0 "$filename";

counter=1;

echo "num_thread,execution-time,iterations" >> "$filename";

export OMP_SCHEDULE=static;

while [ $counter -le 32 ]
do
  export OMP_NUM_THREADS=$counter;
  echo "num threads: $counter";
  time=$(srun --nodes=1 ./a2-omp --n "$1" --m "$2" --epsilon "$epsilon" --max-iterations "$3" --no-verify);
  echo "$time";
  IFS=' ';
  read -ra strarr <<< "$time";
  echo "$counter,${strarr[2]},${strarr[5]}" >> "$filename";
  #echo "${time//$'\t'/','}" >> execution-time-"$1".csv;

  ((counter++));
done
