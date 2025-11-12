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
	if (n == 1) {return ;}
	int m = lastb+1;
	int nlogm = int(n/log2(m));
	if((int(log2(m))%2) == 1)
	{
		nlogm++;
		nlogm = nlogm % n; //this just handles size 2. 
	}

	//step 1 -> create the SRANK arrays. 
	/** We need to make local and global SRANKS.
	 *  Locals should be filled w/ 0s because if not we'll get nulls and be sad
	 * 	Then add in the ranks for n/logm or m/logn
	 *  lastly just gather it all together so everybody has the SRANK arrays
	**/
	int * SRANKA = new int[p*nlogm];
	int * SRANKB = new int[p*nlogm];

	//
	//create the local SRANKS. These will be gathered.
	int * localSRANKA = new int[nlogm](); //() at the end makes it a 0 array.
	int * localSRANKB = new int[nlogm](); 

	//we rank every n/logm or m/logn index of A and B. This call does this. 
	localSRANKA[nlogm-1]=Rank(b, 0,lastb,a[(nlogm*(my_rank+1)-1)%n]);

	localSRANKB[nlogm-1]=Rank(a, 0,lasta,b[((nlogm)*(my_rank+1)-1)%n]);
	//mpi gather the SRANK arrays on all procs.
	MPI_Allgather(localSRANKA, nlogm, MPI_INT, SRANKA, nlogm, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(localSRANKB, nlogm, MPI_INT, SRANKB, nlogm, MPI_INT, MPI_COMM_WORLD);
	delete [] localSRANKA;
	delete [] localSRANKB;

	//GET RID OF BAD PROCS. GO AWAY
	p = log2(n);
	if(my_rank >=p)
	{
		int * bigShape = new int [n+m]();
		MPI_Allreduce(bigShape, output, n+m, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		return;
	}


	//step 2 -> Get the Shapes. Every proc has 2 shapes which we'll call A and B shapes. 
	/**get the A shape for my proc. Each proc will get an A shape and a B shape. 
	*A shape A-side STARTS at a[SRANKB[(m/logn *my_rank)-1]]  || RIGHT
	*A shape A-side ENDS at a[(n/logm)*(my_rank+1)-1]	|| RIGHT
	*A shape B-side STARTS at b[n/logm*my_rank]	 || RIGHT
	*A shape B-side ENDS at b[SRANKA[((n/logm)*(my_rank+1))-1]-1]	|| RIGHT
	*smerge these. 
	**/
	int aShapeAStart;
	if(((nlogm *my_rank)-1) < 0)
	{
		aShapeAStart = 0; //case where we start too low. so we just set to 0. 
	}
	else
	{
		aShapeAStart = SRANKB[int((nlogm *my_rank)-1)];
	}
	int aShapeBStart = nlogm*my_rank;
	int aShapeAEnd = int(((nlogm)*(my_rank+1)-1) - aShapeAStart); //need - start for smerge
	int aShapeBEnd = (SRANKA[int(((nlogm)*(my_rank+1))-1)]-1) - aShapeBStart;
	int shapeASize = aShapeAEnd + aShapeBEnd + 2;
	
	/**get the B shape for my proc.  
	 * B shape B-side STARTS at b[SRANKA[(my_rank+1)*(n/logm)-1]+1]
	 * B shape A-side STARTS at a[(my_rank+1)*(n/logm)]
	 * B shape B-sied ENDS at b[(m/logn)*(my_rank+1)-1]
	 * B shape A-side ENDS at a[SRANKB[(my_rank+1)*(m/logn)-1]-1]
	**/
	
	int bShapeBStart = SRANKA[int((my_rank+1)*(nlogm)-1)];
	int bShapeAStart = (my_rank+1)*(nlogm);
	int bShapeBEnd = ((nlogm)*(my_rank+1)-1)	- bShapeBStart; // - for smerge
	int bShapeAEnd = (SRANKB[int((my_rank+1)*(nlogm)-1)]-1)		-bShapeAStart;

	if(my_rank == 2)
	{
		cout << "TEST FOR RANK 2" << endl;
		cout << aShapeBStart << endl;
		cout << aShapeAStart << endl;
		cout << aShapeBEnd << endl;
		cout << aShapeAEnd << endl;
	}
	
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
	
	
	if((my_rank == p-1) && (p%2==1))
	{
		//weird case where we dont have a b shape because theres an odd num of shapes
		cout << "we made it here " << endl;
		aShapeAEnd = (n-1) - aShapeAStart;
		aShapeBEnd = (n-1) - aShapeBStart;
		smerge(a + aShapeAStart, b + aShapeBStart, aShapeAEnd, aShapeBEnd,bigShape+startIndex);
	}
	else
	{
		//be normal
		//shape a merge
		smerge(a + aShapeAStart, b + aShapeBStart, aShapeAEnd, aShapeBEnd,bigShape+startIndex);
		//shape b merge
		smerge(b + bShapeBStart,a + bShapeAStart,bShapeBEnd,bShapeAEnd,bigShape + shapeASize +startIndex);
	}


	cout << "BIG SHAPE FOR RANK " << my_rank << " - | ";
	for(int i =0; i < n+m; i++)
	{
		cout << bigShape[i] << "| ";
	}
	cout << endl;

	//STEP 4: SHARE IT WE ARE DONE AND NO LONGER SAD 
	MPI_Allreduce(bigShape, output, n+m, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	cout << "|";
	for (int i = 0; i < 4; i++)
	{
		cout << output[i] << "| ";
	}
	cout << endl << endl;
	delete [] bigShape;
	delete [] SRANKA;
	delete [] SRANKB;
}

void mergesort (int * a, int first, int last, int p, int my_rank) //on first pass, last should just be size of array-1. 
{
	if(last-first <7) return;

	if (first >= last) return;

    int mid = (first + last) / 2;
	cout << "run " << last+1 << endl;
    mergesort(a, first, mid, p, my_rank);
	cout << "merged " << first << ", " << mid << endl;
    mergesort(a, mid + 1, last, p, my_rank);
	cout << "merged " << mid+1 << ", " << last << endl;
	
	
	int * tempa = a + first;
	int * tempb = a + mid + 1;
	

	
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
	int size = 32;
	
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
	a1[16] = 4;
	a1[17] = 6;
	a1[18] = 8;
	a1[19] = 9;
	a1[20] = 10;
	a1[21] = 11;
	a1[22] = 15;
	a1[23] = 17;
	a1[24] = 19;
	a1[25] = 20;
	a1[26] = 21;
	a1[27] = 24;
	a1[28] = 25;
	a1[29] = 27;
	a1[30] = 31;
	a1[31] = 32;

	int * b1 = new int [16];
	b1[0] = 1;
	b1[1] = 2;
	b1[2] = 4;
	b1[3] = 5;
	b1[4] = 6;
	b1[5] = 9;
	b1[6] = 12;
	b1[7] = 14;
	b1[8] = 3;
	b1[9] = 7;
	b1[10] = 8;
	b1[11] = 10;
	b1[12] = 11;
	b1[13] = 13;
	b1[14] = 15;
	b1[15] = 16;

	int * c1 = new int [8];
	c1[0] = 1;
	c1[1] = 2;
	c1[2] = 5;
	c1[3] = 6;
	c1[4] = 3;
	c1[5] = 4;
	c1[6] = 7;
	c1[7] = 8;

	if (my_rank == 0)
	{
		for(int i =0; i<16; i++)
		{
			cout << b1[i] << " | ";
		}
		cout << endl;
	}

	mergesort(b1, 0, 15, p, my_rank);

	if(my_rank == 0)
	{
		cout << "FINAL OUTPUT : |";
		for(int i = 0; i<16; i++)
		{
			cout << b1[i] << "| ";
		}
		cout << "END" << endl;
	}

	delete [] a1;
	delete [] b1;

	// Shut down MPI
	MPI_Finalize();

	return 0;
}