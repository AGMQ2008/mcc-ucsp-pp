#include <iostream>
#include <math.h>
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

int main()
{
	// Loading initial data
	initData(N, masses, pos, vel);

	// Clear forces array
	clearForces(N, forces);

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
			cout << "Simulacion en t = " << (delta_t * s) << endl;
			printPosAndVel(N, pos, vel);
		}
		// Reiniciando las fuerzas para la prox simulacion
		for(int q=0; q<N; q++)
		{
			forces[q][X] = forces[q][Y] = 0;
		}

		for(int q=0; q<N; q++)
		{
			// Calcula la fuerza total sobre q
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
		for(int q=0; q<N; q++)
		{
			// Calcula la posicion y velocidad de q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Resultado final\n";
	printPosAndVel(N, pos, vel);

	return 0;
}
