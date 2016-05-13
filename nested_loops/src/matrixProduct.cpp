#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <ctime>

using namespace std;

void matrixInit(int N, int* A)
{
	for (int i = 0; i<N; i++)
	{
		A[i] = rand() % 100;
	}
}

void matrixInit_2D(int N, int** A)
{
	for (int i = 0; i<N; i++)
		for (int j=0; j<N; j++)
		{
			A[i][j] = rand() % 100;
		}
}

void matrixReset(int N, int *A)
{
	for (int i = 0; i<N; i++)
	{
		A[i] = 0;
	}
}

void matrixReset_2D(int N, int **A)
{
	for (int i = 0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			A[i][j] = 0;
		}
}

void matrixMult(int N, int A[], int B[], int* C)
{
	int i, j, k;

	for (i = 0; i<N; i++)
		for (j = 0; j<N; j++)
			for (k = 0; k<N; k++)
				C[i*N + j] += A[i*N + k] * B[k*N + j];
}

void matrixMult_2D_3loops(int N, int **A, int **B, int** C)
{
	int i, j, k;

	for (i = 0; i<N; i++)
		for (j = 0; j<N; j++)
			for (k = 0; k<N; k++)
				C[i][j] += A[i][k] * B[k][j];
}

void matrixMult_2D_6loops(int N, int **A, int **B, int** C)
{
	int numBlocks = 2;
	int sizeBlock = N / numBlocks;

	int **Ctemp;
	Ctemp = new int *[sizeBlock];
	for(int t=0; t<sizeBlock; t++)
		Ctemp[t] = new int[sizeBlock];

	for(int i=0; i<numBlocks; i++)
		for(int j=0; j<numBlocks; j++)
		{
			matrixReset_2D(sizeBlock, Ctemp);

			for (int k = 0; k<sizeBlock; k++)
				for (int ib = 0; ib < sizeBlock; ib++)
					for (int jb = 0; jb<sizeBlock; jb++)
						Ctemp[ib][jb] += A[k][i*sizeBlock+ib] * B[j*sizeBlock+jb][k];

			for (int ib=0; ib<sizeBlock; ib++)
			{
				for(int jb=0; jb<sizeBlock; jb++)
					C[j*sizeBlock+jb][i*sizeBlock+ib] = Ctemp[jb][ib];
			}
		}

}

void matrixPrint(int N, int A[])
{
	int k = N;
	int len = N * N;

	for (int i = 0; i<len; i++)
	{
		if (i == k)
		{
			cout << "\n";
			k += N;
		}
		cout << A[i] << " ";
	}
	cout << endl;
}

void matrixPrint_2D(int N, int **A)
{
	for (int i = 0; i<N; i++)
	{
		for(int j=0; j<N; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}
}

void matrixTest1D()
{
	const int N = 292;		// con un numero mayor stack overflow
	const int len = N * N;

	int A[len], B[len], C[len];

	matrixInit(len, A);
	matrixInit(len, B);
	matrixReset(len, C);

	clock_t begin = clock();
	matrixMult(N, A, B, C);
	clock_t end = clock();

	double elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;

	//cout << "Matriz A: \n";
	//matrixPrint(N, A);
	//cout << endl << "Matriz B: \n";
	//matrixPrint(N, B);
	//cout << endl << "Matriz C: \n";
	//matrixPrint(N, C);

	cout << endl << "Time elapsed of 3 loops: " << elapsed_secs << endl;
}

void matrixTest2D()
{
	const int N = 4;

	int **A, **B, **C;

	A = new int *[N];
	B = new int *[N];
	C = new int *[N];

	for(int i=0; i<N; i++)
	{
		A[i] = new int[N];
		B[i] = new int[N];
		C[i] = new int[N];
	}

	matrixInit_2D(N, A);
	matrixInit_2D(N, B);
	matrixReset_2D(N, C);

	clock_t begin = clock();
	matrixMult_2D_3loops(N, A, B, C);
	clock_t end = clock();

	/*
	cout << "Matriz A: \n";
	matrixPrint_2D(N, A);
	cout << endl << "Matriz B: \n";
	matrixPrint_2D(N, B);
	cout << endl << "Matriz C: \n";
	matrixPrint_2D(N, C);
	*/

	double elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;
	cout << endl << "Time elapsed of 3 loops: " << elapsed_secs << endl;

	matrixPrint_2D(N, C);

	matrixReset_2D(N, C);
	begin = clock();
	matrixMult_2D_6loops(N, A, B, C);
	end = clock();

	elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;
	cout << endl << "Time elapsed of 6 loops: " << elapsed_secs << endl;
	matrixPrint_2D(N, C);
}

int main()
{
	//matrixTest1D();
	matrixTest2D();
	return 0;
}
