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

	cout << "Performing floyds algorithm for matrix of size " << N << endl;
	int k,i,j;
	int ij,ik,kj;

	// Maximum path length is N so we iterate N times
	for(k=0; k<N; k++){
		// Test rows
		for(i=start; i<start+count; i++){
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
					if(data[ij] == 0) data[ij] = FLOYDINF;
					// If our data is smaller, replace it and set the output to be the current path length
					if(data[ik]+data[kj]< data[ij]){
						data[ij] = data[ik]+data[kj];
					}
				}
			}
		}
	}
}

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
	// Broadcast out the matrix width/height
	MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Broadcast out the matrix contents
	MPI_Bcast (data, N*N, MPI_INT, 0, MPI_COMM_WORLD);

	int count = (int) ceil(N/size);
	int start = N;

	FloydsAlgorithm(data,N*N,0,count);
}
// Slave process - receives a request, performs floyd's algorithm, and returns a subset of the data
void Slave(int rank,int S){
	int N;
	MPI_Status status;

	// Receive broadcast of N (the width/height of the matrix)
	MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int num = N;
	int size = N * N;
	int data[size];

	// Receive the matrix
	MPI_Bcast (&data, size, MPI_INT, 0, MPI_COMM_WORLD);

	// Calculate start and count
	int count = (int) ceil(num/S);
	int start = rank * count;
	if((num * start) + (num * count) > size) count = N - start;

	cout << "Process " << rank << "(" << start << "," << count << ")" << endl;

	FloydsAlgorithm(data,num,start,count);
	int index;
	for(int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			index = i*N+j;
			if(data[index] == FLOYDINF)
				cout << 0 << ' ';
			else
				cout << data[index] << ' ';
		}
		cout << endl;
	}

	// Total number of individual items processed
	int total = num * count;
	// Output
	int * out = new int[total];

	// MPI_Send(&out, total, MPI_INT, 0, rank, MPI_COMM_WORLD);
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

	if (rank == 0)
	{  Server(size); }
	else
	{  Slave(rank,size); }
	MPI_Finalize();
}
