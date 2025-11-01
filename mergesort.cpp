#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include <math.h>
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

void pmerge(int * a, int * b, int lasta, int lastb, int p, int my_rank, int * output = NULL)
{
	int n = lasta+1;
	int m = lastb+1;

	//step 1 -> create the SRANK arrays. 
	int * SRANKA = new int[n];
	for(int i = 0; i < n; i++)
	{
		SRANKA[i] = 0;
	}
	int * SRANKB = new int[m];
	for(int i = 0; i < m; i++)
	{
		SRANKB[i] = 0;
	}

	int * localSRANKA = new int[n/p]; //may need to change this. 
	int * localSRANKB = new int[m/p]; 

	//we rank every n/logm or m/logn index of A and B. This call does this. 
	localSRANKA[int(n/log2(m)-1)]=Rank(b, 0,lastb,a[int((n/log2(m))*(my_rank+1)-1)]);

	if (my_rank == 0)
	{
		cout << "the rank - " << Rank(b, 0,lastb,a[int((n/log2(m))*(my_rank+1)-1)]) << endl;
		cout << "AFTER RANK SRANKA: ";
		for(int i =0; i<n/p; i++)
		{
			cout << SRANKA[i] << " | ";
		}
		cout << endl;
	}

	localSRANKB[int((m/log2(n))-1)]=Rank(a, 0,lasta,b[int((m/log2(n))*(my_rank+1)-1)]);
	//mpi gather the SRANK arrays on all procs.
	MPI_Allgather(localSRANKA, n/p, MPI_INT, SRANKA, n/p, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(localSRANKB, m/p, MPI_INT, SRANKB, m/p, MPI_INT, MPI_COMM_WORLD);
	delete [] localSRANKA;
	delete [] localSRANKB;


	//SRANK TEST PRINT
	if(my_rank == 0)
	{
		cout << "SRANKA: |";
		for(int i = 0; i<n; i++)
		{
			cout << SRANKA[i] << " | ";
		}
		cout << endl;

		cout << "SRANKB: |";
		for(int i = 0; i<m; i++)
		{
			cout << SRANKB[i] << " | ";
		}
		cout << endl;
	}
	

	//step 2 -> Get the Shapes. 
	/**get the A shape for my proc. Each proc will get an A shape and a B shape. 
	*A shape A-side STARTS at a[SRANKB[(m/logn *my_rank)-1]]  || RIGHT
	*A shape A-side ENDS at a[(n/logm)*(my_rank+1)-1]	|| RIGHT
	*A shape B-side STARTS at b[n/logm*my_rank]	 || RIGHT
	*A shape B-side ENDS at b[SRANKA[((n/logm)*(my_rank+1))-1]-1]	|| RIGHT
	*smerge these. 
	**/

	/*
	int aShapeAStart;
	if((m/log2(n) *my_rank)-1) < 0)
	{
		aShapeAStart = SRANKB[0];
	}
	else
	{
		aShapeAStart = SRANKB[(m/log2(n) *my_rank)-1];
	}
	int aShapeBStart = n/log2(m)*my_rank;
	int aShapeAEnd = ((n/log2(m))*(my_rank+1)-1)	- aShapeAStart; //need this - for smerge.
	int aShapeBEnd = (SRANKA[((n/log2(m))*(my_rank+1))-1]-1)	- bShapeBStart;
	int ShapeA[aShapeAEnd + aShapeBEnd];
	smerge(a[aShapeAStart], b[aShapeBStart], aShapeAEnd, bShapeBEnd,ShapeA);

	*/

	/**get the B shape for my proc.  
	 * B shape B-side STARTS at b[SRANKA[(my_rank+1)*(n/logm)-1]+1]
	 * B shape A-side STARTS at a[(my_rank+1)*(n/logm)]
	 * B shape B-sied ENDS at b[(m/logn)*(my_rank+1)-1]
	 * B shape A-side ENDS at a[SRANKB[(my_rank+1)*(m/logn)-1]-1]
	**/

	/*
	int bShapeBStart = SRANKA[(my_rank+1)*(n/log2(m))-1]+1;
	int bShapeAStart = (my_rank+1)*(n/log2(m));
	int bShapeBEnd = ((m/log2(n))*(my_rank+1)-1)	- bShapeBStart; // - for smerge
	int bShapeAEnd = (SRANKB[(my_rank+1)*(m/log2(n))-1]-1)		-bShapeAStart;
	int ShapeB[bShapeAEnd + bShapeBEnd];
	smerge(b[bShapeBStart],a[bShapeAStart],bShapeBEnd,bShapeAEnd,ShapeB);
	*/
	delete SRANKA;
	delete SRANKB;
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
	int size = 16;
	
	int * a1 = new int [size];
	a1[0] = 1;
	a1[1] = 2;
	a1[2] = 3;
	a1[3] = 5;
	a1[4] = 7;
	a1[5] = 12;
	a1[6] = 13;
	a1[7] = 14;
	a1[8] = 16;
	a1[9] = 18;
	a1[10] = 22;
	a1[11] = 23;
	a1[12] = 26;
	a1[13] = 28;
	a1[14] = 29;
	a1[15] = 30;
	
	int * b1 = new int[size];
	b1[0] = 4;
	b1[1] = 6;
	b1[2] = 8;
	b1[3] = 9;
	b1[4] = 10;
	b1[5] = 11;
	b1[6] = 15;
	b1[7] = 17;
	b1[8] = 19;
	b1[9] = 20;
	b1[10] = 21;
	b1[11] = 24;
	b1[12] = 25;
	b1[13] = 27;
	b1[14] = 31;
	b1[15] = 32;

	if (my_rank == 0)
	{
		cout << "rank test - " << Rank(b1, 0, 15, 14); 

		for(int i =0; i<size; i++)
		{
			cout << a1[i] << " | ";
		}
		cout << endl;

		for(int i =0; i<size; i++)
		{
			cout << b1[i] << " | ";
		}
		cout << endl;
	}
	

	int * c1 = new int [32];
	pmerge(a1, b1, 15, 15, p,my_rank,c1);
	
	delete [] a1;
	delete [] b1;
	delete [] c1;

	// Shut down MPI
	MPI_Finalize();

	return 0;
}