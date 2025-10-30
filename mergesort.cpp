#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include "mpi.h" // message passing interface           - an api allowing processors to communicate with each other. 
using namespace std;

int Rank(int * a, int first, int last, int valToFind)
{
	// Binary Search
	int rank = 0;
	int low = first;
	int high = last - first;
	while(low <= high)
	{
		int mid = low + (high - low) / 2;
		
		if(a[mid] == valToFind)
		{
			return rank;
		}
		else if(a[mid] < valToFind)
		{
			low = mid + 1;
			rank++;
		}
		else
		{
			high = mid - 1;
			rank++;
		}
	}
	return rank;
}

void smerge(int * a, int * b, int lasta, int lastb, int * output = NULL)
{	
	//the two fingers - each array only increments once it finds a valid answer. 
	//lastb is the last index in b. lasta is the last index in a. 
	int x =0; //index in c
	int i = 0; //index in a
	int j = 0; //index in b
	int segmentSize = lasta + lastb + 2;
	int c [segmentSize]; // size of the chunk of the array. last = last index aka size -1. 
	while(i <= lasta && j <= lastb)
	{
		if(a[i] <= b[j])
		{
			c[x] = a[i];
			i++;
		}
		else
		{
			c[x] = b[j];
			j++;
		}
		x++;
	}
	
	//there might be remaining elements
	while(i <= lasta)
	{
		c[x] = a[i];
		i++;
		x++;
	}
	while(j <= lastb)
	{
		c[x] = b[j];
		j++;
		x++;
	}
	
	//save the sorted array to output. bc this is a pointer, itll just save onto the big array at its chunk position. 
	for(int k = 0; k<segmentSize; k++)
	{
		output[k] = c[k];
	}
}

void mergesort (int * a, int first, int last) //on first pass, last should just be size of array-1. 
{
	if (first >= last){return;} //if we get to 1 left. 
    
	// New
	int lasta = (last - first) / 2 + first;
	int lastb = last;
	int firstb = lasta + 1;
	int firsta = first;
	int mida = (lasta - firsta) / 2 + firsta;
	int midb = (lastb - firstb) / 2 + firstb;
	int sizeA = lasta - firsta;
	int sizeB = lastb - firstb;
	
	// Initializations need changed. Maybe int sizeA = firsta - lasta + 1. Same for sizeB
	int * tempa = new int[sizeA];
	int * tempb = new int[sizeB];

	mergesort(tempa, firsta, mida - 1);
	mergesort(tempa, mida, lasta);
	smerge(tempa + firsta, tempa + mida, mida - firsta - 1, last - mida, tempa + firsta);

	mergesort(tempb, firstb, midb - 1);
	mergesort(tempb, midb, lastb);
	smerge(tempb + firstb, tempb + midb, midb - firstb - 1, last - midb, tempb + firstb);

	int * c = new int[last + 1];
	// pmerge next
}

// New compile and run commands for MPI!
// mpicxx -o blah file.cpp
// mpirun -q -np 32 blah

// Linux commands to help you manage processes
// ps
// kill -9

int main (int argc, char * argv[]) {

	int my_rank;			// my CPU number for this process
	int p;					// number of CPUs that we have
	int source;				// rank of the sender
	int dest;				// rank of destination
	int tag = 0;			// message number
	char message[100];		// message itself
	MPI_Status status;		// return status for receive
	
	// Start MPI
	MPI_Init(&argc, &argv);
	
	// Find out my rank!
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	// Find out the number of processes!
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	// THE REAL PROGRAM IS HERE
    //get user number for array size
	int arrSize;
	cout << "Please enter the size of the array: " << endl;
	cin >> arrSize;
	
	int a[arrSize];
	srand(time(nullptr));
	for(int i = 0; i<arrSize; i++)
	{
		a[i] = rand()%1000 +1;
	}
	
	//debug - display the initial array
	for(int i =0; i<arrSize; i++)
	{
		cout << a[i] << " | ";
	}
	cout << endl;
	
	//actually sort it. 
	mergesort(a,0,arrSize-1);
	
	//display the sorted array. 
	for(int i =0; i<arrSize; i++)
	{
		cout << a[i] << " | ";
	}
	cout << endl;
	
	// Shut down MPI
	MPI_Finalize();

	return 0;
}