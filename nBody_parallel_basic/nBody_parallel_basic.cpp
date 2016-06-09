#include <iostream>
#include <math.h>
#include <omp.h>
#include "funcBasic.h"

using namespace std;

// Global variables
#define N 			100	// Numero de particulas
#define delta_t		1	// Intervalo de tiempo
#define numSim		3	// Numero de simulaciones
#define tipoRes		0	// 0: Resultado final, 1: Resultado por iteracion

#define G 	6.673

double masses[N];
double pos[N][2], vel[N][2];
double forces[N][2];

int main()
{
	// Loading initial data
	initData(N, masses, pos, vel);
	//loadDataFromFile("particles.txt", N, masses, pos, vel);

	// Clear forces array
	clearForces(N, forces);

	double x_diff, y_diff, dist, dist_cubed;
	int numThreads;

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
		numThreads = omp_get_num_threads();
#pragma omp for schedule(static, N / numThreads)
		for(int q=0; q<N; q++)
		{
			forces[q][X] = forces[q][Y] = 0;

			// Calculate total force in q
			for(int k=0; k<N; k++)
			{
				if(k != q)
				{
					x_diff = pos[q][X] - pos[k][X];
					y_diff = pos[q][Y] - pos[k][Y];
					dist = sqrt(x_diff*x_diff + y_diff*y_diff);
					dist_cubed = dist * dist * dist;
					forces[q][X] -= G * masses[q]*masses[k] / dist_cubed * x_diff;
					forces[q][Y] -= G * masses[q]*masses[k] / dist_cubed * y_diff;
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

	return 0;
}
