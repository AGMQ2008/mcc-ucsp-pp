#include <iostream>
#include <math.h>
#include <pthread.h>
#include "funcBasic.h"

using namespace std;

// Global variables
#define N 			10000	// Numero de particulas
#define delta_t		1	// Intervalo de tiempo
#define numSim		3	// Numero de simulaciones
#define tipoRes		0	// 0: Resultado final, 1: Resultado por iteracion

#define G 	6.673 //E-11

#define numThreads 4

double masses[N];
double pos[N][2], vel[N][2];
double forces[N][2];

typedef struct arg_struct {
	int rank;
	int ini;
	int fin;
	int inc;
};

int counter = 0;
pthread_mutex_t mutex;
pthread_cond_t cond_var;

void *calculateForces(void *arguments)
{
	int ini, fin, inc;
	arg_struct *args = (arg_struct *) arguments;

	ini = args->ini;
	fin = args->fin;
	inc = args->inc;

	double x_diff, y_diff, dist, dist_cubed;

	for(int q=ini; q<fin; q+=inc)
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
	pthread_mutex_lock(&mutex);
	counter++;
	if(counter == numThreads)
	{
		counter = 0;
		pthread_cond_broadcast(&cond_var);
	}
	else
	{
		//cout << "Hilo "<< args->rank << " pasa a estado de espera..\n";
		while(pthread_cond_wait(&cond_var, &mutex) != 0) ;
	}
	pthread_mutex_unlock(&mutex);

	//cout << "Actualizacion de pos y vel por hilo "<< args->rank << "\n";
	for(int q=ini; q<fin; q+=inc)
	{
		// Calcule position and velocity of q
		pos[q][X] += delta_t * vel[q][X];
		pos[q][Y] += delta_t * vel[q][Y];
		vel[q][X] += delta_t / masses[q] * forces[q][X];
		vel[q][Y] += delta_t / masses[q] * forces[q][Y];
	}

	return NULL;
}

arg_struct createBlock(int threadId, int blockSize, int numParticles)
{
	arg_struct data;
	data.rank = threadId;
	data.ini = threadId * blockSize;
	data.fin = min(threadId * blockSize + blockSize, numParticles);
	data.inc = 1;
	return data;
}

void Loop_schedule(int numThr, int nPar)
{
	pthread_t threads[numThr];
	int block = nPar / numThr;

	arg_struct args[numThr];

	for(int i=0; i<numThr; i++)
	{
		args[i] = createBlock(i, block, nPar);
		pthread_create(&threads[i], NULL, calculateForces, (void *)&args[i]);
	}

	for(int i=0; i<numThr; i++)
		pthread_join(threads[i], NULL);
}

int main()
{
	// Loading initial data
	initData(N, masses, pos, vel);

	// Clear forces array
	clearForces(N, forces);

	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
			cout << "Simulation in t = " << (delta_t * s) << endl;
			printPosAndVel(N, pos, vel);
		}
		counter = 0;
		Loop_schedule(numThreads, N);
	}
	cout << "Final result\n";
	//printPosAndVel(N, pos, vel);

	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&cond_var);

	return 0;
}
