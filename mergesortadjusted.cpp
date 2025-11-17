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

    if (a[last] < valToFind) {
        return last + 1;}
	return low;
}

void pmerge(int * a, int * b, int lasta, int lastb, int p, int my_rank, int * output = NULL)
{
	int n = lasta+1;
    int m = lastb+1;
	if (n == 1 || m == 1) 
    {
        smerge(a, b, lasta, lastb, output);
        return ;
    }
	int nlogm = int(n/int(log2(m)));
    int mlogn = int(m/int(log2(n)));
    int numShapes = (int(n/nlogm) + int(m/mlogn));
    if( numShapes % 2 == 1)
    {
        //numShapes--;
    }

	//step 1 -> create the SRANK arrays. 
	/** We need to make local and global SRANKS.
	 *  Locals should be filled w/ 0s because if not we'll get nulls and be sad
	 * 	Then add in the ranks for n/logm or m/logn
	 *  lastly just gather it all together so everybody has the SRANK arrays
	**/
	int * SRANKA = new int[n];
	int * SRANKB = new int[m];

	//
	//create the local SRANKS. These will be gathered.
	int * localSRANKA = new int[n](); //() at the end makes it a 0 array.
	int * localSRANKB = new int[m](); 

	//we rank every n/logm or m/logn index of A and B. This call does this. 
    for(int i =my_rank; i < (n/nlogm); i+=p)
    {
        localSRANKA[(nlogm*(i+1)-1)]=Rank(b, 0,lastb,a[(nlogm*(i+1)-1)%n]);
    }

    for(int i =my_rank; i < (m/mlogn); i+=p)
    {
        localSRANKB[((mlogn)*(i+1)-1)]=Rank(a, 0,lasta,b[((mlogn)*(i+1)-1)%m]);
    }

	//mpi gather the SRANK arrays on all procs.
    MPI_Allreduce(localSRANKA, SRANKA, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(localSRANKB, SRANKB, m, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	delete [] localSRANKA;
	delete [] localSRANKB;

    // we need to parse both A and B and find where segments stop and start. 
    //For odd shapes, the A side stops at n/logm. for even shapes, the A side stops at SRANKB((logn/m)-1)
    //for odd shapes, the B side stops at SRANKA[(logn/m) -1]. For even shapes, the B side stops at logn/m.

    //we need a way to set this up so that it works on any number of procs. This means if there are 24 procs and 8 shapes,
    //we only need 8. It also means that if there are less than 8 procs, we need to have the procs handle every pth.
    // we could set this up with a loop that goes until we are over the amount of shapes (determined above) and then increments by i+p.
    // once in the loop, find our given shapes, as well as where they start in the main overall arr. 
    // then go ahead and copy them to bigShape using smerge. Each arr should have it's own local bigShape instantiated to 0.
    //one thing that may make this easier is to share 2 global arrays for the stop points of A and B. Each proc finds one. We can use this later. 

    int * aPoints = new int [numShapes+1]();
    int * bPoints = new int [numShapes+1]();


    for(int i = my_rank; i< int(m/mlogn); i+=p)
    {
        //evens
        int aStopPoint = SRANKB[(mlogn * (i+1)) -1];
        aPoints[i+1+int(n/nlogm)] = aStopPoint;
        int bStopPoint = mlogn * (i+1);
        //deal with a odd number of total shapes
        if(i == (int(m/mlogn))-1 && bStopPoint != m)
        {
            bStopPoint = m;
        }
        bPoints[i+1] = bStopPoint;
    }

    for (int i = my_rank; i < int (n/nlogm); i+=p)
    {
        //odds 
        int aStopPoint = nlogm * (i+1);
        //deal with a odd number of total shapes
        if(i == (int(n/nlogm))-1 && aStopPoint != n)
        {
            aStopPoint = n;
        }
        aPoints[i+1] = aStopPoint;
        int bStopPoint = SRANKA[(nlogm * (i+1)) -1];
        bPoints[i+1+int(m/mlogn)] = bStopPoint;
    }

    //cout << "POINTS FOR RANK " << my_rank +1 << " : A - |";
    //cout << "numshapes+1 is " << numShapes+1 << endl;
    //gather the stop points together. We can just do allreduce and sum. 
    int * tempA = new int [numShapes+1];
    int * tempB = new int [numShapes+1];
    MPI_Allreduce(aPoints, tempA, numShapes+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(bPoints, tempB, numShapes+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    delete [] aPoints;
    delete [] bPoints;


    //now we just need to smerge each half of the globalPoints arrays.
    int * globalAPoints = new int [numShapes+1];
    globalAPoints[0] = 0;
    int * globalBPoints = new int [numShapes+1];
    globalBPoints[0] = 0;

    smerge(tempA+1, tempA+(int(n/nlogm))+1, int(n/nlogm)-1, int(m/mlogn)-1, globalAPoints+1);
    smerge(tempB+1, tempB+(int(m/mlogn))+1, int(m/mlogn)-1, int(n/nlogm)-1, globalBPoints+1);
    
    if(my_rank == 0)
	{
		cout << "APOINTS : |";
        for (int i = 0; i < numShapes+1; i++)
        {
            cout << globalAPoints[i] << " | ";
        }
        cout << endl;

        cout << "BPOINTS : |";
        for (int i = 0; i < numShapes+1; i++)
        {
            cout << globalBPoints[i] << " | ";
        }
        cout << endl;
    }

    //then do same process but make the shapes using the start and stop points. 
    int * bigShape = new int [n+m](); //() makes it 0s.
    for(int i = my_rank; i < numShapes; i+=p)
    {
        // our starting point for each shape is going to be the combination of the prev stop points
        int startPoint = globalAPoints[i] + globalBPoints[i];
        if(i == 0){startPoint = 0;}
        int aStart = globalAPoints[i];
        int bStart = globalBPoints[i];
        int aEnd = (globalAPoints[i+1]-1) - globalAPoints[i];
        int bEnd = (globalBPoints[i+1]-1) - globalBPoints[i];

        smerge(a + aStart, b + bStart, aEnd, bEnd,bigShape + startPoint);
    }
    
    //then lastly combine down like before. 
    MPI_Allreduce(bigShape, output, n+m, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    cout << "FAILTEST: " << a[0] << " , " << b[0] << ", " << lasta << ", " << lastb << endl;
    delete [] globalAPoints;
    delete [] globalBPoints;
    delete [] bigShape;
	delete [] SRANKA;
	delete [] SRANKB;
}

void mergesort (int * a, int first, int last, int p, int my_rank) //on first pass, last should just be size of array-1. 
{
	//if(last-first <1) return;

	if (first >= last) return;

    int mid = (first + last) / 2;
    mergesort(a, first, mid, p, my_rank);
    mergesort(a, mid + 1, last, p, my_rank);
	
	int * tempa = a + first;
	int * tempb = a + mid + 1;
	
	pmerge(tempa, tempb, mid - first, last - mid - 1, p, my_rank, a + first);

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
	a1[0] = 4;
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
	a1[16] = 1;
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
    if (my_rank == 0)
    {
        cout << Rank(a1, 0, 15, 31) << endl;
    }

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

	int * c1 = new int [16];
	c1[0] = 2;
	c1[1] = 3;
	c1[2] = 10;
	c1[3] = 15;
	c1[4] = 5;
	c1[5] = 6;
	c1[6] = 9;
	c1[7] = 11;
	c1[8] = 1;
	c1[9] = 7;
	c1[10] = 12;
	c1[11] = 13;
	c1[12] = 4;
	c1[13] = 8;
	c1[14] = 14;
	c1[15] = 16;

	int * d1 = new int [32];

	mergesort(a1, 0, 30, p, my_rank);
	//pmerge (a1, a1+16, 15, 13, p, my_rank, d1);

	if(my_rank == 0)
	{
		cout << "FINAL OUTPUT : |";
		for(int i = 0; i<30; i++)
		{
			cout << a1[i] << "| ";
		}
		cout << "END" << endl;
	}


	delete [] a1;
	delete [] b1;
	delete [] c1;
    delete [] d1;

	// Shut down MPI
	MPI_Finalize();

	return 0;
}