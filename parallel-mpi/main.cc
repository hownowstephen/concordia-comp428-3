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
#include <fstream>
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
 * @param int start The start of the block
 * @param int end The end of the block
 */
void FloydsAlgorithm(int rank, int *data, int N, int start, int count){

	int k,i,j,k_here;
	int ij,ik;

	int rowk[N];

	for (k=0;k<N;k++) {

		int owner = (int) ceil(k/count);
		// Check if k is owned by this process, if so, calculate the row
		if (rank == owner) {
			k_here = k;// - start;
			for(j=0;j<N;j++)
				rowk[j]=data[k_here*N + j];
		}

		//MPI_Bcast(&k, 1, MPI_INT, 0, MPI_Comm_World);
		MPI_Bcast(rowk,N,MPI_INT,owner,MPI_COMM_WORLD);

		for(i=start;i<start+count;i++){
			for(j=0;j<N;j++){

				ij = i * N + j;
				ik = i * N + k;

				if(i == j){
					data[ij] = 0;
				}else{
					if(data[ij] == 0) data[ij] = SMINF;
					data[ij]= min(data[ij], data[ik]+rowk[j]);
				}
			}
		}
	}
}

void Server(int size,char * file){
	// Perform dispatch of all requests
	// Need to: Broadcast data, send each process a start/count pair for their requirements
	MPI_Status status;

	FILE *I_in;
	// Load in the Adjacency matrix to test
	ifstream M_in(file, ios::in);
	int N,tmp,index;
	M_in >> N;

	// Generate the dataset
	int data[N*N];
	for (int y = 0; y < N; y++)
		for (int x = 0; x < N; x++){
			   M_in >> tmp;
			   data[y*N + x] = tmp;
		}
	// Finally, print the result
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
	// Broadcast out the matrix width/height
	MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Broadcast out the matrix contents
	MPI_Bcast (data, N*N, MPI_INT, 0, MPI_COMM_WORLD);

	int count = (int) ceil(N/size);
	int start = N;
	int total = N*N;

	FloydsAlgorithm(0,data,N,0,count);

	int t[total];
	for(int p=1;p<size;p++){
		MPI_Recv(&t, total, MPI_INT, p, 0, MPI_COMM_WORLD,&status);
		for(int v=0;v<total;v++){
			data[v] = max(data[v],t[v]);
		}
	}
	// Finally, print the result
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

	// Perform my transformations
	FloydsAlgorithm(rank,data,num,start,count);

	// Send my data
	MPI_Send(data,size,MPI_INT,0,0,MPI_COMM_WORLD);
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

	char * file;

	// Take a filename as a param
	if(argc > 1){
		 file = argv[1];
	}else{
		cout << "Please supply a filename" << endl;
		MPI_Finalize();
		return 1;
	}

	if (rank == 0)
	{
		double startTime = getClock();
		Server(size,file);
		cout << "Time: "<< getClock() - startTime << endl;
	}
	else
	{  Slave(rank,size); }
	MPI_Finalize();
}
