#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include "mpi.h" // message passing interface           - an api allowing processors to communicate with each other. 
using namespace std;

int rank(int * a, int first, int last, int valToFind)
{
	//needs to be Binary Search -> goes to logn time instead of n time. 
	int rank = 0;
	for(int i = first; i<= last; i++)
	{
		if(valToFind > a[i])
		{
			rank++;
		}
		else
		{
			return rank;
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
        
	int mid = first + ((last - first+1) / 2);//gets the middle point (aka first pos in RHS)
	mergesort(a, first, mid-1);
	mergesort(a, mid, last);
	smerge(a + first, a + mid, mid-first-1, last-mid, a+first); //a start, b start, a last, b last. 

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