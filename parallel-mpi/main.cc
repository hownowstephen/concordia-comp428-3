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

// Define a smaller infinity constant
#define SMINF 999999

char name[MPI_MAX_PROCESSOR_NAME];


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
 * @param int rank The rank of the current process
 * @param int* data A square adjacency matrix
 * @param int N Size of a single dimension of the matrix
 * @param int i The row to test
 */
void FloydsAlgorithm(int rank, int *data, int N, int start, int count){

	cout << "Performing floyds algorithm for matrix of size " << N << endl;
	cout << "Start " << start << " count " << count << endl;
	int k,i,j,k_here;
	int ij,ik;

	int rowk[N];

	for (k=0;k<N;k++) {
		// Check if k is owned by this process, if so, calculate the row
		if (k >= start && k < start+count) {
			k_here = k - start;
			for(j=0;j<N;j++)
				rowk[j]=data[k_here*N + j];
		}

		//MPI_Bcast(&k, 1, MPI_INT, 0, MPI_Comm_World);
		MPI_Bcast(rowk,N,MPI_INT,rank,MPI_COMM_WORLD);

		for(i=start;i<start+count;i++){
			for(j=0;j<N;j++){

				ij = i * N + j;
				ik = i * N + k;

				if(data[ij] == 0) data[ij] = SMINF;
				if(i == j) data[ij] = 0;
				data[ij]= min(data[ij], data[ik]+rowk[j]);
			}
		}
	}

	if(start == 0){
		int index;
		for(int i=0;i<N;i++){
			for (int j=0;j<N;j++){
				index = i*N+j;
				if(data[index] == SMINF)
					cout << 0 << ' ';
				else
					cout << data[index] << ' ';
			}
			cout << endl;
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
					1, 0, 1, 0,
				  };
	// Broadcast out the matrix width/height
	MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Broadcast out the matrix contents
	MPI_Bcast (data, N*N, MPI_INT, 0, MPI_COMM_WORLD);

	int count = (int) ceil(N/size);
	int start = N;

	FloydsAlgorithm(0,data,N,0,count);
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

	FloydsAlgorithm(rank,data,num,start,count);

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
