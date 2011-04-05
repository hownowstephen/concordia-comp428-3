/**
 * Parallel implementation of Floyds all-pairs shortest path algorithm
 * using the HP-MPI C Library
 *
 * Resources used in development:
 * http://www.mcs.anl.gov/~itf/dbpp/text/node35.html
 * http://www.scribd.com/doc/23340655/11/Floyd%E2%80%99s-sequential-algorithm
 *
 *  @author Stephen Young (st_youn@encs.concordia.ca)
 *  @course Comp428/4
 *
 */

#include <iostream>
#include <algorithm>
#include <mpi.h>
#include <sys/time.h>
using namespace std;

char name[MPI_MAX_PROCESSOR_NAME];

// Define a constant for infinity
#define INFINITY 99999

/**
 * GetClock
 * Used for benchmarking performance
 */
double getClock()
{
    struct timeval tp;
    struct timezone tzp;

    gettimeofday(&tp, &tzp);

    return (double)tp.tv_sec+((double)tp.tv_usec/1000000.0);
}


/**
 * FloydsAlgorithm (parallelized)
 * Performs a test of an adjacency matrix (2-dimensional represented as single dimension)
 * to find the shortest distance between two nodes on the directed graph
 * @param int* data A square adjacency matrix
 * @param int N Size of a single dimension of the matrix
 * @param int i The row to test
 */
void FloydsAlgorithm(int *data, int N, int i){

	int k,j;
	int ij,ik,kj;

	// Output will be same size as data
	int out[N];

	// Maximum path length is N so we iterate N times
	for(k=0; k<N; k++){
		// Test columns
		for (j=0; j<N; j++){
			// Resolve some indices
			ij = i * N + j;
			ik = i * N + k;
			kj = k * N + j;

			// If i == j, the nodes are the same, so the distance is zero
			if(i==j){
				data[ij] = 0;
			}else{
				// Use an arbitrarily large number to test against
				if(data[ij] == 0) data[ij] = INFINITY;
				// If our data is smaller, replace it and set the output to be the current path length
				if(data[ik]+data[kj]< data[ij]){
					data[ij] = data[ik]+data[kj];
					// Only need index j
					out[j] = k;
				}
			}
		}
	}
}

// Stub for server process
void Server(int size, int rank){
	// Perform dispatch of all requests
	// Need to: Broadcast data, send each process a start/count pair for their requirements
	Slave(rank); // Server partakes as a slave as well
}
// Stub for slave process
void Slave(int rank){
	// Handles all received orders
}

/**
 * Main function
 * Initializes the required mpi communication layer
 * and dispatches both the server and the slave processes
 * The server will also act as a slave to ensure that all the processors are
 * busy.
 */
int main(void){
	int size, rank, len;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(name, &len);

	if (rank == 0)
	{  Server(size,rank); }
	else
	{  Slave(rank); }
	MPI_Finalize();
}
