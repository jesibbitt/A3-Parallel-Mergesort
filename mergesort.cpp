#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include <math.h>
#include "mpi.h" // message passing interface           - an api allowing processors to communicate with each other. 
using namespace std;

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

int Rank(int * a, int first, int last, int valToFind)
{
	// Binary Search
	int low = first;
	int high = last;
	while(low < high)
	{
		int mid = low + (high - low) / 2;
		
		if(a[mid] < valToFind)
		{
			low = mid + 1;
		}
		else
		{
			high = mid;
		}
	}
	return low;
}

void pmerge(int * a, int * b, int lasta, int lastb, int p, int my_rank, int * output = NULL)
{
	int n = lasta+1;
	int m = lastb+1;

	//step 1 -> create the SRANK arrays. 
	/** We need to make local and global SRANKS.
	 *  Locals should be filled w/ 0s because if not we'll get nulls and be sad
	 * 	Then add in the ranks for n/logm or m/logn
	 *  lastly just gather it all together so everybody has the SRANK arrays
	**/
	int * SRANKA = new int[n];
	int * SRANKB = new int[m];

	//create the local SRANKS. These will be gathered.
	int * localSRANKA = new int[n/p](); //() at the end makes it a 0 array.
	int * localSRANKB = new int[m/p](); 

	//we rank every n/logm or m/logn index of A and B. This call does this. 
	localSRANKA[int(n/log2(m)-1)]=Rank(b, 0,lastb,a[int((n/log2(m))*(my_rank+1)-1)]);

	localSRANKB[int((m/log2(n))-1)]=Rank(a, 0,lasta,b[int((m/log2(n))*(my_rank+1)-1)]);
	
	//mpi gather the SRANK arrays on all procs.
	MPI_Allgather(localSRANKA, n/p, MPI_INT, SRANKA, n/p, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(localSRANKB, m/p, MPI_INT, SRANKB, m/p, MPI_INT, MPI_COMM_WORLD);
	delete [] localSRANKA;
	delete [] localSRANKB;

	//step 2 -> Get the Shapes. Every proc has 2 shapes which we'll call A and B shapes. 

	/**get the A shape for my proc. Each proc will get an A shape and a B shape. 
	*A shape A-side STARTS at a[SRANKB[(m/logn *my_rank)-1]]  || RIGHT
	*A shape A-side ENDS at a[(n/logm)*(my_rank+1)-1]	|| RIGHT
	*A shape B-side STARTS at b[n/logm*my_rank]	 || RIGHT
	*A shape B-side ENDS at b[SRANKA[((n/logm)*(my_rank+1))-1]-1]	|| RIGHT
	*smerge these. 
	**/

	int aShapeAStart;
	if(((m/log2(n) *my_rank)-1) < 0)
	{
		aShapeAStart = SRANKB[0]; //case where we start too low. so we just set to 0. 
	}
	else
	{
		aShapeAStart = SRANKB[int((m/log2(n) *my_rank)-1)];
	}
	int aShapeBStart = n/log2(m)*my_rank;
	int aShapeAEnd = int(((n/log2(m))*(my_rank+1)-1) - aShapeAStart); //need - start for smerge
	int aShapeBEnd = (SRANKA[int(((n/log2(m))*(my_rank+1))-1)]-1) - aShapeBStart;
	int shapeASize = aShapeAEnd + aShapeBEnd + 2;
	
	/**get the B shape for my proc.  
	 * B shape B-side STARTS at b[SRANKA[(my_rank+1)*(n/logm)-1]+1]
	 * B shape A-side STARTS at a[(my_rank+1)*(n/logm)]
	 * B shape B-sied ENDS at b[(m/logn)*(my_rank+1)-1]
	 * B shape A-side ENDS at a[SRANKB[(my_rank+1)*(m/logn)-1]-1]
	**/
	
	int bShapeBStart = SRANKA[int((my_rank+1)*(n/log2(m))-1)];
	int bShapeAStart = (my_rank+1)*(n/log2(m));
	int bShapeBEnd = ((m/log2(n))*(my_rank+1)-1)	- bShapeBStart; // - for smerge
	int bShapeAEnd = (SRANKB[int((my_rank+1)*(m/log2(n))-1)]-1)		-bShapeAStart;

	//STEP 3: MAKE OUR COMBINED SHAPE FOR OUR PROC
	/** We need to make a big shape of A+B so that we can send both back
	 * We do this by using a bigShape array, and calling smerge so that it returns-
	 * -each shape to its appropriate spot in the big shape- 
	 * -(B after A, at the index it will be in for the final sol)
	 * This is literally just 2 smerge calls.
	**/

	//these are special cases where one side of the B shape has nothing.
	if (bShapeAEnd < 0) // if nothing in the a side. 
	{
		bShapeAStart = 0; //this makes sure we dont search for something that doesnt exist. 
		bShapeAEnd = -1; //THIS -1 is VERY VERY important for smerge. 
	}
	else if (bShapeBEnd < 0) // if nothing is in the b side
	{
		bShapeBStart = 0;
		bShapeBEnd = -1;
	}

	//make the big shape
	int * bigShape = new int [n+m]();
	int startIndex = min(a[aShapeAStart]-1, b[aShapeBStart]-1);

	//shape a merge
	smerge(a + aShapeAStart, b + aShapeBStart, aShapeAEnd, aShapeBEnd,bigShape+startIndex);
	//shape b merge
	smerge(b + bShapeBStart,a + bShapeAStart,bShapeBEnd,bShapeAEnd,bigShape + shapeASize +startIndex);

	//STEP 4: SHARE IT WE ARE DONE AND NO LONGER SAD 
	MPI_Allreduce(bigShape, output, n+m, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	delete [] bigShape;
	delete [] SRANKA;
	delete [] SRANKB;
}

void mergesort (int * a, int first, int last) //on first pass, last should just be size of array-1. 
{
	if (first >= last) return;

    int mid = (first + last) / 2;
    mergesort(a, first, mid, p, my_rank);
    mergesort(a, mid + 1, last, p, my_rank);
	
	
	int * tempa = a + first;
	int * tempb = a + mid + 1;
	
	cout << "pre-pmerge" << endl;
	
	pmerge(tempa, tempb, mid - first, last - mid - 1, p, my_rank, a + first);
	
	cout << "pmerge done" << endl;

	//smerge(a + first, a + mid + 1, mid - first, last - mid - 1, a + first);
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
		cout << "rank test - " << Rank(b1, 0, 15, 14) << endl; 

		for(int i =0; i<size; i++)
		{
			cout << a1[i] << " | ";
		}
		cout << endl;

	if (my_rank == 0)
	{
		cout << "rank test - " << Rank(b1, 0, 15, 14) << endl; 

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

	if(my_rank == 0)
	{
		cout << "FINAL OUTPUT : |";
		for(int i = 0; i<32; i++)
		{
			cout << c1[i] << "| ";
		}
		cout << endl;
	}
		for(int i =0; i<size; i++)
		{
			cout << b1[i] << " | ";
		}
		cout << endl;
	}
	

	int * c1 = new int [32];
	pmerge(a1, b1, 15, 15, p,my_rank,c1);

	if(my_rank == 0)
	{
		cout << "FINAL OUTPUT : |";
		for(int i = 0; i<32; i++)
		{
			cout << c1[i] << "| ";
		}
		cout << endl;
	}

	
	
	delete [] a1;
	delete [] b1;
	delete [] c1;

	// Shut down MPI
	MPI_Finalize();

	return 0;
}