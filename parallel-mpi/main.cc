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
#include <math.h>
#include <sys/time.h>
using namespace std;

char name[MPI_MAX_PROCESSOR_NAME];

// Define a constant for infinity
#define FLOYDINF 999999

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
void FloydsAlgorithm(int *data, int N, int start, int count){

	int k,j,i;
	int ij,ik,kj;

	// Output will be same size as data
	int out[N];

	// Maximum path length is N so we iterate N times
	for(k=0; k<N; k++){
		// Test columns
		for (j=0; j<N; j++){
			// Process the set of columns indicated
			for(i=start;i<start+count;i++){
				// Resolve some indices
				ij = i * N + j;
				ik = i * N + k;
				kj = k * N + j;

				// If i == j, the nodes are the same, so the distance is zero
				if(i==j){
					data[ij] = 0;
				}else{
					// Use an arbitrarily large number to test against
					if(data[ij] == 0) data[ij] = FLOYDINF;
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
}

// Stub for server process
void Server(int size){
	// Perform dispatch of all requests
	// Need to: Broadcast data, send each process a start/count pair for their requirements

	int N = 4;
	/** @var sample A sample adjacency matrix */
	int data[16]= { 0, 1, 0, 0,
					0, 0, 1, 1,
					0, 0, 0, 0,
					1, 0, 1, 0
				  };

	cout << "Initializing Floyd's Algorithm";

	// Broadcast out the matrix width/height
	MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Broadcast out the matrix contents
	MPI_Bcast (data, N*N, MPI_INT, 1, MPI_COMM_WORLD);

	int count = ceil(N/size);
	int start = N;

	// Send directives to each processor of what to process
	for (int dest = 1; dest < size; ++dest){
		MPI_Send (&start, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
		MPI_Send (&count, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
		start += count;
		// Handle the last processor potentially having less rows
		if(start + count > N) count = N - start;

	}

	FloydsAlgorithm(data,N,0,ceil(N/size));
}
// Slave process - receives a request, performs floyd's algorithm, and returns a subset of the data
void Slave(int rank){
	int N,start,count;
	MPI_Status status;

	// Receive broadcast of N (the width/height of the matrix)
	MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int * data = new int[N*N];
	// Receive the matrix
	MPI_Bcast (&data, N*N, MPI_INT, 1, MPI_COMM_WORLD);

	// Receive directives for processing
	MPI_Recv (&start, 1, MPI_INT, rank, 0, MPI_COMM_WORLD,&status);
	MPI_Recv (&count, 1, MPI_INT, rank, 1, MPI_COMM_WORLD,&status);

	FloydsAlgorithm(data,N,start,count);

	// Total number of individual items processed
	int total = N * count;
	// Output
	int * out = new int[total];
	// Populate the output
	int c = 0;
	for(int i=start;i<start+count;i++){
		for(int j=0;j<N;j++){
			out[c] = data[i*j];
			cout << out[c] << " ";
		}
		cout << endl;
	}
	MPI_Send(&out, total, MPI_INT, SERVER, rank, MPI_COMM_WORLD);
}

/**
 * Main function
 * Initializes the required mpi communication layer
 * and dispatches both the server and the slave processes
 * The server will also act as a slave to ensure that all the processors are
 * busy.
 */
int main(int argc, char * argv[]){
	int size, rank, len;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(name, &len);

	cout << "Rank " << rank << " checking in" << endl;

	if (rank == 0)
	{  Server(size); }
	else
	{  Slave(rank); }
	MPI_Finalize();
}
