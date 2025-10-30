#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include "mpi.h" // message passing interface           - an api allowing processors to communicate with each other. 
using namespace std;

void pmerge(int * a, int * b, int lasta, int lastb, int * output = NULL)
{
	int n = lasta+1;
	int m - lastb+1;

	//step 1 -> create the SRANK arrays. 
	int SRANKA[n];
	int SRANKB[m];
	int localSRANKA [n/p]; //may need to change this. 
	int localSRANKB [m/p]; 

	//we rank every n/logm or m/logn index of A and B. This call does this. 
	localSRANKA[(n/log2(m))-1]=Rank(b, 0,lastb,a[(n/log2(m))*(my_rank+1)-1]);
	localSRANKB[(m/log2(n))-1]=Rank(a, 0,lasta,b[(m/log2(n))*(my_rank+1)-1]);
	//mpi gather the SRANK arrays on all procs.
	MPI_Allgather(localSRANKA, n/p, MPI_INT, SRANKA, n/p, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(localSRANKB, m/p, MPI_INT, SRANKB, m/p, MPI_INT, MPI_COMM_WORLD);

	//SRANK TEST PRINT
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
}

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
    
	// New
	int lasta = (last - first) / 2 + first;
	int lastb = last;
	int firstb = lasta + 1;
	int firsta = first;
	int mida = (lasta - firsta) / 2 + firsta;
	int midb = (lastb - firstb) / 2 + firstb;
	
	// Initializations need changed. Maybe int sizeA = firsta - lasta + 1. Same for sizeB
	int * tempa = a[first];
	int * tempb = a[firstb];

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