#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "mpi/mpi.h"

#include "a2-helpers.hpp"

// tags that indicate if we send data upwards or downwards
#define UPWARDS_TAG 1
#define DOWNWARDS_TAG 2

#define WORLD MPI_COMM_WORLD

using namespace std;

int main(int argc, char **argv)
{
    int max_iterations = 1000;
    double epsilon = 1.0e-3;
    bool verify = true, print_config = false;

    // default values for M rows and N columns
    int N = 12;
    int M = 12;

    process_input(argc, argv, N, M, max_iterations, epsilon, verify, print_config);

    int numprocs, rank;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(WORLD, &numprocs);
    MPI_Comm_rank(WORLD, &rank);

    if ( print_config )
        std::cout << "Configuration: m: " << M << ", n: " << N << ", max-iterations: " << max_iterations << ", epsilon: " << epsilon << std::endl;

    auto time_mpi_1 = MPI_Wtime();

    // START: the main part of the code that needs to use MPI/OpenMP timing routines 
    //  MPI hint: remember to initialize MPI first 
    
    int i, j;
    double diffnorm;
    int iteration_count = 0;

    // change M to local size for horizontal fragmentation
    int old_M = M;
    int local_M = M / numprocs + 2; // +2 for padding rows

    // from and to are the indices for each data chunk
    int from = rank * M;
    // make last process handle remaining rows
    int to = (rank == numprocs - 1) ? old_M - 1 : from + (M - 1);

    Mat U(local_M, N); // for MPI: use local sizes with MPI, e.g., recalculate M and N
    Mat W(local_M, N); // for MPI: use local sizes with MPI, e.g., recalculate M and N

    // Init & Boundary
    for (i = 0; i < local_M; ++i) {
        for (j = 0; j < N; ++j) {
            W[i][j] = U[i][j] = 0.0;
        }

        W[i][0] = U[i][0] = 0.05; // left side
        W[i][N-1] = U[i][N-1] = 0.1; // right side
    }

    for (j = 0; j < N; ++j) {
        W[1][j] = U[1][j] = 0.02; // top
        W[local_M - 2][j] = U[local_M - 2][j] = 0.2; // bottom
    }
    // End init

    MPI_Request request_upwards, request_downwards;
    MPI_Status  status_upwards, status_downwards;

    iteration_count = 0;
    do
    {
        iteration_count++;
        diffnorm = 0.0;

        // Compute new values (but not on boundary) 
        for (i = 1; i < local_M - 1; ++i)
        {
            for (j = 1; j < N - 1; ++j)
            {
                W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
            }
        }

        // Only transfer the interior points
        for (i = 1; i < local_M - 1; ++i)
            for (j = 1; j < N - 1; ++j)
                U[i][j] = W[i][j];

        diffnorm = sqrt(diffnorm); // all processes need to know when to stop

        if (rank == 0) {
            //std::cout << "Iteration: " << iteration_count << ", diffnorm: " << diffnorm << ", Epsilon: " << epsilon << std::endl;
        }

        if (epsilon >= diffnorm) {
            cout << "Rank: " << rank << ", Iteration: " << iteration_count << ", diffnorm: " << diffnorm << ", Epsilon: " << epsilon << endl;
        }

        if (rank != (numprocs - 1)) {
            // receive from below
            // send downwards

            MPI_Irecv(U[local_M-1], N, MPI_DOUBLE, rank + 1, UPWARDS_TAG, WORLD, &request_upwards);
            MPI_Send(U[local_M-2], N, MPI_DOUBLE, rank + 1, DOWNWARDS_TAG, WORLD);
            MPI_Wait(&request_upwards, &status_upwards);

            //MPI_Isend(U[local_M-2], N, MPI_DOUBLE, rank + 1, DOWNWARDS_TAG, MPI_COMM_WORLD, &request_downwards);
            //MPI_Recv(U[local_M-1], N, MPI_DOUBLE, rank + 1, UPWARDS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // send bottom row downwards for all chunks except the last one
        if (rank != 0) {
            // receive from above
            // send upwards

            MPI_Irecv(U[0], N, MPI_DOUBLE, rank - 1, DOWNWARDS_TAG, WORLD, &request_downwards);
            MPI_Send(U[1], N, MPI_DOUBLE, rank - 1, UPWARDS_TAG, WORLD);
            MPI_Wait(&request_downwards, &status_downwards);

            //MPI_Isend(U[1], N, MPI_DOUBLE, rank - 1, UPWARDS_TAG, MPI_COMM_WORLD, &request_upwards);
            //MPI_Irecv(U[0], N, MPI_DOUBLE, rank - 1, DOWNWARDS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } while (epsilon <= diffnorm && iteration_count < max_iterations);

    Mat final_local_U(local_M - 2, N);

    // copy local U matrix to final_local_U matrix without boundary rows
    for (i = 0; i < local_M - 2; ++i) {
        for (j = 0; j < N; ++j) {
            final_local_U[i][j] = U[i + 1][j];
        }
    }

    // final matrix
    Mat bigU(M, N);

    // gather all local parts of the final_local_U matrix for rank 0
    MPI_Gather(&final_local_U(0, 0), (local_M - 2) * N, MPI_DOUBLE, &bigU(0, 0), (local_M - 2) * N, MPI_DOUBLE, 0, WORLD);

    if (rank != 0){
        MPI_Finalize();
        return 0;
    }

    auto time_mpi_2 = MPI_Wtime(); // change to MPI_Wtime() / omp_get_wtime()

    // Print time measurements 
    cout << "Elapsed time: "; 
    cout << std::fixed << std::setprecision(4) << chrono::duration<double>(time_mpi_2 - time_mpi_1).count(); // remove for MPI/OpenMP
    // cout << std::fixed << std::setprecision(4) << (time_2 - time_1); // modify accordingly for MPI/OpenMP
    cout << " seconds, iterations: " << iteration_count << endl; 
 
    // verification     
    if ( verify ) {
        Mat U_sequential(M, N); // init another matrix for the verification

        int iteration_count_seq = 0;
        heat2d_sequential(U_sequential, max_iterations, epsilon, iteration_count_seq); 

        // Here we need both results - from the sequential (U_sequential) and also from the OpenMP/MPI version, then we compare them with the compare(...) function 
        cout << "Verification: " << ( U.compare(U_sequential) && iteration_count == iteration_count_seq ? "OK" : "NOT OK") << std::endl;
    }

    MPI_Finalize();

    return 0;
}