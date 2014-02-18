#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

////////////////////////////////////////////////////////////
//
// This is my implementation of sequence alignment for class
// COT 5405, spring 2013. Anyone can use this code for any
// purpose for free. However, there is no guarantee for the  
// quality or reliability of this code. 
//
// If you find any bug, you can let me know, but I don't have
// any plan to fix it. However, if you fixed the bug and 
// mailed your copy to me, I would like to update this code.
// 
// Bo Yang(bonny95@gmail.com), March 2013.
//
///////////////////////////////////////////////////////////

using namespace std;

const int MATCH_COST=0;
const int MISMATCH_COST=2;
const int GAP_COST=1;

int memcnt=0; // counting memory allocated

struct BTNode
{
	int x; // x axis
	int y; // y axis
	char type; // 'm'-match,'s'-mismatch,'g'-gap
	char nt; // valid neucleotide: x-sequence x, y-sequence y, b-both 
};

BTNode* P;
int pidx=0; // index of P

struct Region // region in a sequence
{
	int start; 
	int end;
};

//
// Generate random sequence
// 
char* GenRandSeq(char* seq, int n)
{
	string nt="ATCG";
	for(int i=0;i<n;++i)
	{
		int idx=rand()%4;
		seq[i]=nt[idx];
	}

	return seq;
}

int min(int a,int b,int c)
{
	return (a<b?a:b)<c?(a<b?a:b):c;
}

void Merge(BTNode* input, long p, long r)
{
	long mid = floor((p + r) / 2);
	long i1 = 0;
	long i2 = p;
	long i3 = mid + 1;

	// Temp array
	BTNode* temp=new BTNode[r-p+1];
	memcnt+=sizeof(BTNode)*(r-p+1);

	// Merge in sorted form the 2 arrays
	while ( i2 <= mid && i3 <= r )
		if ( input[i2].x < input[i3].x )
			temp[i1++] = input[i2++];
		else
			temp[i1++] = input[i3++];

	// Merge the remaining elements in left array
	while ( i2 <= mid )
		temp[i1++] = input[i2++];

	// Merge the remaining elements in right array
	while ( i3 <= r )
		temp[i1++] = input[i3++];

	// Move from temp array to master array
	for ( int i = p; i <= r; i++ )
		input[i] = temp[i-p];

	delete [] temp;
}

// Merge sort the BTNode array by the first index(element x)
//
// inputs:
//	p - the start index of array input
//	r - the end index of array input
void Merge_sort(BTNode* input, long p, long r)
{
	if ( p < r )
	{
		long mid = floor((p + r) / 2);
		Merge_sort(input, p, mid);
		Merge_sort(input, mid + 1, r);
		Merge(input, p, r);
	}
}

//
// Insert new node into P
//
void InsertNodetoP(int x, int y)
{
	bool isok=true;
	for(int i=0;i<pidx;++i)
	{
		if(x==P[i].x && y==P[i].y)
			isok=false;
	}
	if(isok)
	{
		P[pidx].x=x;
		P[pidx].y=y;
		pidx++;
	}
}

//
// Align seqences
//
// Recurrence:
// 	OPT(i,j)=min{   
// 			alpha(i,j)+OPT(i-1,j-1), 
// 			delta+OPT(i-1,j), 
// 			delta+OPT(i,j-1)
// 		}
// 	where mismatch cost alpha=2, match cost alpha=0, and gap cost delta=1.
//
// Parameters:
//	X,Y - sequences
//	m - length of sequence X
//	n - length of sequence Y
//
int Alignment(char* X, Region& rx, char* Y, Region& ry)
{
	int m=rx.end-rx.start+1;
	int n=ry.end-ry.start+1;

	int** S=new int*[m+1]; // temp array to store scores of alignment
	for(int i=0;i<m+1;++i)
		S[i]=new int[n+1];
	memcnt+=sizeof(int)*(m+1)*(n+1);

	for(int i=0;i<=m;++i)
		S[i][0]=i*GAP_COST;
	for(int j=0;j<=n;++j)
		S[0][j]=j*GAP_COST;
	
	for(int j=1;j<=n;++j)
	{
		for(int i=1;i<=m;++i)
		{
			int alpha;
			if(X[i-1]==Y[j-1])
				alpha=MATCH_COST;
			else
				alpha=MISMATCH_COST;

			S[i][j]=min(alpha+S[i-1][j-1],GAP_COST+S[i-1][j],GAP_COST+S[i][j-1]);
		}
	} // end of for

	// Trace back 
	int ix=m;
	int iy=n;
	int bt_len=m+n;// length of array to record tracing back
	int cnt=bt_len-1; 	
	BTNode* bt=new BTNode[bt_len]; // record backtrace path
	memcnt+=sizeof(BTNode)*bt_len;
	while(ix>=1 && iy>=1)
	{
		bt[cnt].x=ix-1;
		bt[cnt].y=iy-1;
		InsertNodetoP(rx.start+ix-1,ry.start+iy-1);
				
		if(S[ix][iy]==S[ix-1][iy-1]+MATCH_COST && X[ix-1]==Y[iy-1]) // match
		{
			bt[cnt].type='m';
			bt[cnt].nt='b';
			ix--;
			iy--;
		} else if(S[ix][iy]==S[ix-1][iy-1]+MISMATCH_COST) { // mismatch
			bt[cnt].type='s';
			bt[cnt].nt='b';
			ix--;
			iy--;
		} else if(S[ix][iy]==S[ix-1][iy]+GAP_COST) { // gap
			bt[cnt].type='g';
			bt[cnt].nt='x';
			ix--;
		} else if(S[ix][iy]==S[ix][iy-1]+GAP_COST) { // gap
			bt[cnt].type='g';
			bt[cnt].nt='y';
			iy--;
		} 	
		cnt--;
	}
	while(iy>0)
	{
		bt[cnt].x=ix-1;
		bt[cnt].y=iy-1;
		InsertNodetoP(rx.start+ix-1,ry.start+iy-1);
	
		if(S[ix][iy]==S[ix][iy-1]+GAP_COST) { // gap
			bt[cnt].type='g';
			bt[cnt].nt='y';
			iy--;
		}
		cnt--;
	}
	while(ix>0)
	{
		bt[cnt].x=ix-1;
		bt[cnt].y=iy-1;
		InsertNodetoP(rx.start+ix-1,ry.start+iy-1);
		
		if(S[ix][iy]==S[ix-1][iy]+GAP_COST) { // gap
			bt[cnt].type='g';
			bt[cnt].nt='x';
			ix--;
		}
		cnt--;
	}

/*	// Print the alignment
  	cnt++;
	// TEST ONLY
	for(int i=cnt;i<bt_len;++i)
	{
		cout<<"("<<bt[i].x<<","<<bt[i].y<<") ";
	}
	cout<<endl;

	// Print sequence X
	for(int i=cnt;i<bt_len;++i)
	{
		if(bt[i].nt=='x'||bt[i].nt=='b')
		{
			if(bt[i].x<0)
				cout<<'-';
			else
				cout<<X[bt[i].x];
		} else {
			cout<<'-';
		}
	}
	cout<<endl;
	// print middle line
	for(int i=cnt;i<bt_len;++i)
	{
		if(bt[i].type=='m')
			cout<<'|';
		else
			cout<<' ';
	}
	cout<<endl;
	// print sequence Y
	for(int i=cnt;i<bt_len;++i)
	{
		if(bt[i].nt=='y'||bt[i].nt=='b')
		{
			if(bt[i].y<0)
				cout<<'-';
			else
				cout<<Y[bt[i].y];
		} else {
			cout<<'-';
		}
	}
	cout<<endl; */

	delete [] bt;

	int score=S[m][n];

	for(int i=0;i<m+1;++i)
		delete [] S[i];
	delete [] S;

	return score;
}

// 
// Space efficient alignment, which caculate the length of the 
// shortest path from (0,0) to (i,j) and only can returns the optimal value.
//
int SpaceEfficientAlignment(int** S, char* X, int m, char* Y, int n)
{
	for(int i=0;i<=m;++i)
	{
		S[i][0]=i*GAP_COST;
	}

	for(int j=1;j<=n;++j)
	{
		S[0][1]=j*GAP_COST;
		for(int i=1;i<=m;++i)
		{
			int alpha;
			if(X[i-1]==Y[j-1])
				alpha=MATCH_COST;
			else
				alpha=MISMATCH_COST;

			S[i][1]=min(alpha+S[i-1][0],GAP_COST+S[i-1][1],GAP_COST+S[i][0]);
		}

		// move column 1 of S to column 0 to make room for next iteration
		for(int i=0;i<=m;++i)
			S[i][0]=S[i][1];
	}

	return S[m][1];
}

// 
// Backward space efficient alignment, which caculate the length of the 
// shortest path from (i,j) to (m,n) and only can returns the optimal value.
//
int BackwardSpaceEfficientAlignment(int** S, char* X, int m, char* Y, int n)
{
	for(int i=m;i>=0;--i)
	{
		S[i][1]=(m-i)*GAP_COST;
	}

	for(int j=n-1;j>=0;--j)
	{
		S[m][0]=(n-j)*GAP_COST;
		for(int i=m-1;i>=0;--i)
		{
			int alpha;
			//if(X[m-i-1]==Y[j])
			if(X[i]==Y[j])
				alpha=MATCH_COST;
			else
				alpha=MISMATCH_COST;

			S[i][0]=min(alpha+S[i+1][1],GAP_COST+S[i+1][0],GAP_COST+S[i][1]);
		}

		// move column 0 of S to column 1 to make room for next iteration
		for(int i=0;i<=m;++i)
			S[i][1]=S[i][0];
	}

	return S[0][0];
}

// 
// Find out the q that minimize the score f(q,n/2)+g(q,n/2)
//
// Inputs:
// 	sf - score vector of forward space efficient alignment
// 	sb - score vector of backward space efficient alignment
// 	m - length of the score vector
//
int FindMinScore(int** sf,int** sb,int m)
{
	int min=sf[1][1]+sb[1][0];
	int ret=1;
	for(int q=2;q<=m;++q)
	{
		if(min>sf[q][1]+sb[q][0])
		{
			min=sf[q][1]+sb[q][0];
			ret=q-1; // ret should be the position in the sequence, so minus 1
		}
	}

	return ret;
}

//
// Sequence alignemnt using divide and conquer
//
void DivideConquerAlignment(char* X, Region& rx, char* Y, Region& ry)
{
	int m=rx.end-rx.start+1;
	int n=ry.end-ry.start+1;

	if(m<=2 || n<=2)
	{
		Alignment(&X[rx.start],rx,&Y[ry.start],ry);
//		cout<<"Call Alignment"<<endl; // TEST ONLY
		return;
	}
	
	int** score_forward=new int*[m+1]; // (m+1)x2 matrix to store scores of alignment
	int** score_backward=new int*[m+1]; // (m+1)x2 matrix to store scores of alignment
	for(int i=0;i<m+1;++i)
	{
		score_forward[i]=new int[2];
		score_backward[i]=new int[2];
	}
	memcnt+=sizeof(int)*2*(m+1)*2;

	// Get the middle points in path
	Region tempry;
	tempry.start=ry.start;
	tempry.end=(ry.start+ry.end)/2;
	SpaceEfficientAlignment(score_forward,&X[rx.start],m,&Y[tempry.start],tempry.end-tempry.start+1);
	tempry.start=(ry.start+ry.end)/2+1;
	tempry.end=ry.end;
	BackwardSpaceEfficientAlignment(score_backward,&X[rx.start],m,&Y[tempry.start],tempry.end-tempry.start+1);

	int q=FindMinScore(score_forward,score_backward,m)+rx.start;
//	cout<<"q="<<q<<", n/2="<<(ry.start+ry.end)/2<<endl; // TEST ONLY
	InsertNodetoP(q,(ry.start+ry.end)/2);

	// Divide and Conquer
	Region temprx;
	temprx.start=rx.start;
	temprx.end=q;
	tempry.start=ry.start;
	tempry.end=(ry.start+ry.end)/2;
	if(temprx.start<=temprx.end && tempry.start<=tempry.end)
		DivideConquerAlignment(X,temprx,Y,tempry);
	temprx.start=q+1;
	temprx.end=rx.end;
	tempry.start=(ry.start+ry.end)/2+1;
	tempry.end=ry.end;
	if(temprx.start<=temprx.end && tempry.start<=tempry.end)
		DivideConquerAlignment(X,temprx,Y,tempry);

	for(int i=0;i<m+1;++i)
	{
		delete [] score_forward[i];
		delete [] score_backward[i];
	}
	delete [] score_forward;
	delete [] score_backward;
}

// 
// Backtrace and print the alignment
//
void Backtrace(char* X,char* Y)
{
	Merge_sort(P,0,pidx-1); // Sort nodes in path according to element x
	// Adjust all nodes in path by element y
	int prev_x=-1,prev_y=-1;
	BTNode temp;
	for(int i=0;i<pidx;++i)
	{
		if(prev_y>P[i].y)
		{
			temp=P[i-1];
			P[i-1]=P[i];
			P[i]=temp;
		}
		prev_y=P[i].y;
	}

//	for(int i=0;i<pidx;++i)
//		cout<<"("<<P[i].x<<","<<P[i].y<<")"<<" "; //TEST ONLY
//	cout<<endl; // TEST ONLY

	// Print sequence X
	prev_x=-1;
	for(int i=0;i<pidx;++i)
	{
		if(P[i].x<0 || P[i].x==prev_x)
			cout<<'-';
		else
			cout<<X[P[i].x];
		prev_x=P[i].x;
	}
	cout<<endl;

	// Print middle line
	prev_x=-1;
	prev_y=-1;
	for(int i=0;i<pidx;++i)
	{
		if(X[P[i].x]==Y[P[i].y] && prev_x!=P[i].x && prev_y!=P[i].y)
			cout<<'|';
		else
			cout<<' ';
		prev_x=P[i].x;
		prev_y=P[i].y;
	}
	cout<<endl;

	// Print sequence Y
	prev_y=-1;
	for(int i=0;i<pidx;++i)
	{
		if(P[i].y<0 || P[i].y==prev_y)
			cout<<'-';
		else
			cout<<Y[P[i].y];
		prev_y=P[i].y;
	}
	cout<<endl;

}

//
// MAIN START HERE
//
int main(int argc, char** argv)
{
	// initialize random seed
	srand(time(NULL));

	char* output=NULL;	// output file
	int seqlen=0;		// sequence length
	int c;

	while ((c = getopt (argc, argv, "n:o:")) != -1)
	{
		switch (c)
           	{
		case 'n':
             		seqlen = atoi(optarg); 
             		break;
		case 'o':
             		output = optarg; 
             		break;
           	case '?':
             		if (optopt == 'n')
              			fprintf (stderr, "Option -%c requires an argument.\n\n", optopt);
             		else if (isprint (optopt))
               			fprintf (stderr, "Unknown option `-%c'.\n\n", optopt);
             		else
               			fprintf (stderr,
                        		"Unknown option character `\\x%x'.\n\n",
                        		optopt);
             		return 1;
           	default:
             		abort ();
           	}
	}

	int m=seqlen;
	int n=seqlen-rand()%(seqlen/5);

	char* ref=new char[m];
	char* seq=new char[n];
	memcnt+=sizeof(char)*(m+n);
	
	// generate random sequences
	ref=GenRandSeq(ref,m);
	seq=GenRandSeq(seq,n);

	P=new BTNode[m+n];
	memcnt+=sizeof(BTNode)*(m+n);
	
	cout<<"seq 1:"<<ref<<endl; // TEST ONLY
	cout<<"seq 2:"<<seq<<endl; // TEST ONLY

	Region rx, ry;
	rx.start=0;
	rx.end=m-1;
	ry.start=0;
	ry.end=n-1;
//	Alignment(ref,rx,seq,ry); // TEST ONLY

	DivideConquerAlignment(ref,rx,seq,ry);

	Backtrace(ref,seq);

	delete [] ref;
	delete [] seq;
	delete [] P;

	// Record the space costs
	fstream fs;
	fs.open("memory_cost", fstream::out|fstream::app);
	fs<<m<<"bp: "<<memcnt<<endl;
	fs.close();

	return 0;
}
