#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <ctime>

using namespace std;

void matrixFill_2D(int N, int** A)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = 1;
}

void matrixZeros_2D(int N, int **A)
{
	for (int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			A[i][j] = 0;
}

void matrixPrint_2D(int N, int **A)
{
	for (int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}
}

void matrixMult_2D_3loops(int N, int **A, int **B, int** C)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			for (int k = 0; k < N; k++)
				C[i][j] += A[i][k] * B[k][j];
}

void matrixMult_2D_3loopsC(int N, int **A, int **B, int** C)
{
	for (int i = 0; i < N; i++)
		for (int k = 0; k < N; k++)
			for (int j = 0; j < N; j++)
				C[i][j] += A[i][k] * B[k][j];
}

int BLOCK_SIZE = 16;

void matrixMult_2D_5loops(int N, int **A, int **B, int** C)
{
	for (int i0 = 0; i0 < N; i0 += BLOCK_SIZE)
		for(int j0 = 0; j0 < N; j0 += BLOCK_SIZE)
			for(int i = i0; i < min(i0 + BLOCK_SIZE, N); i++)
				for(int j = j0; j < min(j0 + BLOCK_SIZE, N); j++)
					for(int k = 0; k < N; k++)
						C[i][j] += A[i][k] * B[k][j];
}

void matrixMult_2D_5loopsC(int N, int **A, int **B, int** C)
{
	for (int i0 = 0; i0 < N; i0 += BLOCK_SIZE)
		for(int j0 = 0; j0 < N; j0 += BLOCK_SIZE)
			for(int i = i0; i < min(i0 + BLOCK_SIZE, N); i++)
				for(int k = 0; k < N; k++)
					for(int j = j0; j < min(j0 + BLOCK_SIZE, N); j++)
						C[i][j] += A[i][k] * B[k][j];
}

void matrixMult_2D_6loops(int N, int **A, int **B, int** C)
{
	for (int i = 0; i < N; i += BLOCK_SIZE)
		for (int k = 0; k < N; k += BLOCK_SIZE)
			for (int j = 0; j < N; j += BLOCK_SIZE)
				for (int iInner = i; iInner < min(i + BLOCK_SIZE, N); iInner++)
					for (int kInner = k; kInner < min(k + BLOCK_SIZE, N); kInner++)
						for (int jInner = j; jInner < min(j + BLOCK_SIZE, N); jInner++)
							C[iInner][jInner] += A[iInner][kInner] * B[kInner][jInner];
}

void matrixTest2D()
{
	const int N = 1000;

	int **A, **B, **C;

	clock_t begin, end;
	double elapsed_secs;

	A = new int *[N];
	B = new int *[N];
	C = new int *[N];

	for(int i = 0; i < N; i++)
	{
		A[i] = new int[N];
		B[i] = new int[N];
		C[i] = new int[N];
	}

	matrixFill_2D(N, A);
	matrixFill_2D(N, B);
	matrixZeros_2D(N, C);

	begin = clock();
	matrixMult_2D_3loops(N, A, B, C);
	end = clock();

	elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;
	cout << "Matrix " << N << " x " << N <<" -- Time elapsed of 3 loops normal: " << elapsed_secs << " secs" << endl;
	//matrixPrint_2D(N, C);

	matrixZeros_2D(N, C);
	begin = clock();
	matrixMult_2D_3loopsC(N, A, B, C);
	end = clock();

	elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;
	cout << "Matrix " << N << " x " << N <<" -- Time elapsed of 3 loops by column: " << elapsed_secs << " secs" << endl;
	//matrixPrint_2D(N, C);

	matrixZeros_2D(N, C);
	begin = clock();
	matrixMult_2D_5loops(N, A, B, C);
	end = clock();

	elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;
	cout << "Matrix " << N << " x " << N <<" -- Time elapsed of 5 loops normal: " << elapsed_secs << " secs" << endl;
	//matrixPrint_2D(N, C);

	matrixZeros_2D(N, C);
	begin = clock();
	matrixMult_2D_5loopsC(N, A, B, C);
	end = clock();

	elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;
	cout << "Matrix " << N << " x " << N <<" -- Time elapsed of 5 loops by column: " << elapsed_secs << " secs" << endl;
	//matrixPrint_2D(N, C);

	matrixZeros_2D(N, C);
	begin = clock();
	matrixMult_2D_6loops(N, A, B, C);
	end = clock();

	elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;
	cout << "Matrix " << N << " x " << N <<" -- Time elapsed of 6 loops: " << elapsed_secs << " secs" << endl;
	//matrixPrint_2D(N, C);
}

int main()
{
	matrixTest2D();
	return 0;
}
