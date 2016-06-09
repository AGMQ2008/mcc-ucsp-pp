#include <iostream>
#include <math.h>
#include "omp.h"
#include "funcBasic.h"

using namespace std;

// Global variables
#define N 			4	// Numero de particulas
#define delta_t		1	// Intervalo de tiempo
#define numSim		3	// Numero de simulaciones
#define tipoRes		1	// 0: Resultado final, 1: Resultado por iteracion

#define G 	6.673

double masses[N];
double pos[N][2], vel[N][2];
double forces[N][2];

void nBody_parallel_reduced_locks()
{
	double forces[N][2];

	omp_lock_t *locks;
	locks = new omp_lock_t[N];

	for(int i = 0; i < N; i++)
	{
		forces[i][X] = forces[i][Y] = 0;
		omp_init_lock(&locks[i]);
	}

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

#pragma omp parallel
	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
#pragma omp single
			{
			cout << "Simulation in t = " << (delta_t * s) << endl;
			printPosAndVel(N, pos, vel);
			}
		}
		// Reiniciando las fuerzas para la prox simulacion
#pragma omp for
		for(int q=0; q<N; q++)
		{
			forces[q][X] = forces[q][Y] = 0.0;
		}
#pragma omp for
		for(int q=0; q<N; q++)
		{
			// Calculate total force in q
			for(int k=q+1; k<N; k++)
			{
				x_diff = pos[q][X] - pos[k][X];
				y_diff = pos[q][Y] - pos[k][Y];
				dist = sqrt(x_diff*x_diff + y_diff*y_diff);
				dist_cubed = dist * dist * dist;
				force_qkX = G * masses[q]*masses[k] / dist_cubed * x_diff;
				force_qkY = G * masses[q]*masses[k] / dist_cubed * y_diff;

				omp_set_lock(&locks[q]);
				forces[q][X] += force_qkX;
				forces[q][Y] += force_qkY;
				omp_unset_lock(&locks[q]);

				omp_set_lock(&locks[k]);
				forces[k][X] -= force_qkX;
				forces[k][Y] -= force_qkY;
				omp_unset_lock(&locks[k]);
			}
		}
#pragma omp for
		for(int q=0; q<N; q++)
		{
			// Calcule position and velocity of q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Final result\n";
	printPosAndVel(N, pos, vel);

	for(int i=0; i<N; i++)
		omp_destroy_lock(&locks[i]);
}

void nBody_parallel_reduced_simple()
{
	double forces[N][DIM];

	for(int i = 0; i < N; i++)
	{
		forces[i][X] = forces[i][Y] = 0;
	}

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

#pragma omp parallel
	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
#pragma omp single
			{
			cout << "Simulation in t = " << (delta_t * s) << endl;
			printPosAndVel(N, pos, vel);
			}
		}
		// Reiniciando las fuerzas para la prox simulacion
#pragma omp for
		for(int q=0; q<N; q++)
		{
			forces[q][X] = forces[q][Y] = 0.0;
		}
#pragma omp for
		for(int q=0; q<N; q++)
		{
			// Calculate total force in q
			for(int k=q+1; k<N; k++)
			{
				x_diff = pos[q][X] - pos[k][X];
				y_diff = pos[q][Y] - pos[k][Y];
				dist = sqrt(x_diff*x_diff + y_diff*y_diff);
				dist_cubed = dist * dist * dist;
				force_qkX = G * masses[q]*masses[k] / dist_cubed * x_diff;
				force_qkY = G * masses[q]*masses[k] / dist_cubed * y_diff;

				forces[q][X] += force_qkX;
				forces[q][Y] += force_qkY;
				forces[k][X] -= force_qkX;
				forces[k][Y] -= force_qkY;
			}
		}
#pragma omp for
		for(int q=0; q<N; q++)
		{
			// Calcule position and velocity of q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Final result\n";
	printPosAndVel(N, pos, vel);
}

void nBody_parallel_reduced_critical()
{
	double forces[N][DIM];

	for(int i = 0; i < N; i++)
	{
		forces[i][X] = forces[i][Y] = 0;
	}

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

#pragma omp parallel
	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
#pragma omp single
			{
			cout << "Simulation in t = " << (delta_t * s) << endl;
			printPosAndVel(N, pos, vel);
			}
		}
		// Reiniciando las fuerzas para la prox simulacion
#pragma omp for
		for(int q=0; q<N; q++)
		{
			forces[q][X] = forces[q][Y] = 0.0;
		}
#pragma omp for
		for(int q=0; q<N; q++)
		{
			// Calculate total force in q
			for(int k=q+1; k<N; k++)
			{
				x_diff = pos[q][X] - pos[k][X];
				y_diff = pos[q][Y] - pos[k][Y];
				dist = sqrt(x_diff*x_diff + y_diff*y_diff);
				dist_cubed = dist * dist * dist;
				force_qkX = G * masses[q]*masses[k] / dist_cubed * x_diff;
				force_qkY = G * masses[q]*masses[k] / dist_cubed * y_diff;

#pragma omp critical
				{
					forces[q][X] += force_qkX;
					forces[q][Y] += force_qkY;
					forces[k][X] -= force_qkX;
					forces[k][Y] -= force_qkY;
				}
			}
		}
#pragma omp for
		for(int q=0; q<N; q++)
		{
			// Calcule position and velocity of q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Final result\n";
	printPosAndVel(N, pos, vel);
}

void nBody_parallel_reduced_locforces()
{
	double forces[N][DIM];

	for(int i = 0; i < N; i++)
	{
		forces[i][X] = forces[i][Y] = 0;
	}

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

	int myRank, numThreads;
	double*** loc_forces;

#pragma omp parallel
	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
#pragma omp single
			{
				cout << "Simulation in t = " << (delta_t * s) << endl;
				printPosAndVel(N, pos, vel);
			}
		}
#pragma omp single
		{
		numThreads = omp_get_num_threads();
		loc_forces = new double**[numThreads];

		for(int i=0; i<numThreads; i++)
		{
			loc_forces[i] = new double*[2];
			for(int j=0; j<N; j++)
			{
				loc_forces[i][j] = new double[2];
				loc_forces[i][j][X] = 0;
				loc_forces[i][j][Y] = 0;
			}
		}
		}
#pragma omp for
		for(int q = 0; q < N; q++)
		{
			forces[q][X] = forces[q][Y] = 0;

			// Calculate total force in q
			for(int k = q+1; k < N; k++)
			{
				x_diff = pos[q][X] - pos[k][X];
				y_diff = pos[q][Y] - pos[k][Y];
				dist = sqrt(x_diff*x_diff + y_diff*y_diff);
				dist_cubed = dist * dist * dist;
				force_qkX = G * masses[q]*masses[k] / dist_cubed * x_diff;
				force_qkY = G * masses[q]*masses[k] / dist_cubed * y_diff;

				myRank = omp_get_thread_num();
				loc_forces[myRank][q][X] += force_qkX;
				loc_forces[myRank][q][Y] += force_qkY;
				loc_forces[myRank][k][X] -= force_qkX;
				loc_forces[myRank][k][Y] -= force_qkY;
			}
		}
#pragma omp for
		for(int q = 0; q < N; q++)
		{
			for(int thread = 0; thread < numThreads; thread++)
			{
				forces[q][X] += loc_forces[thread][q][X];
				forces[q][Y] += loc_forces[thread][q][Y];
			}
		}
#pragma omp for
		for(int q = 0; q < N; q++)
		{
			// Calcule position and velocity of q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Final result\n";
	printPosAndVel(N, pos, vel);
}

int main()
{
	// Loading initial data
	initData(N, masses, pos, vel);

	//nBody_parallel_reduced_simple();
	//nBody_parallel_reduced_critical();
	//nBody_parallel_reduced_locks();
	nBody_parallel_reduced_locforces();
}
