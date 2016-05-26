#include <iostream>
#include <fstream>
#include <string>
#include "omp.h"

using namespace std;


double const G = 6.673;
int const X = 0;
int const Y = 1;


void printPosAndVel(int N, double** pos, double** vel)
{
	cout << "Particle \t\t Position \t\t Velocity\n";
	for(int i=0; i<N; i++)
	{
		cout << (i+1) << "\t\t(" << pos[i][X] << "," << pos[i][Y] << ")\t\t "
			<< "(" << vel[i][X] << "," << vel[i][Y] << ")\n";
	}
}

bool loadDataFromFile(string fileName, int N, double *mass, double** pos, double** vel)
{
	fstream myFile(fileName, std::ios_base::in);
	double temp = 0.0;

	if(myFile.is_open())
	{
		for(int i = 0; i < N; i++)
		{
			// Reading mass
			myFile >> temp;
			mass[i] = temp;

			// Reading position
			myFile >> temp;
			pos[i][X] = temp;
			myFile >> temp;
			pos[i][Y] = temp;

			// Reading velocity
			myFile >> temp;
			vel[i][X] = temp;
			myFile >> temp;
			vel[i][Y] = temp;
		}
		return true;
	}
	else
		cout << "Archivo no encontrado\n";
	return false;
}

void nBody_parallel_reduced_locks(string fileName, int N, int delta_t, int numSim, int tipoRes)
{
	// Creando el vector de masas
	double *masses;
	masses = new double[N];

	// Creando los vectores para almacenar las posiciones y las velocidades
	double **pos, **vel, **forces;
	pos = new double *[N];
	vel = new double *[N];
	forces = new double*[N];

	omp_lock_t *locks;
	locks = new omp_lock_t[N];

	for(int i = 0; i < N; i++)
	{
		pos[i] = new double[2];
		vel[i] = new double[2];
		forces[i] = new double[2];
		forces[i][X] = forces[i][Y] = 0;

		omp_init_lock(&locks[i]);
	}

	if(!loadDataFromFile(fileName, N, masses, pos, vel))
		return;

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

#pragma omp parallel
	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
#pragma omp single
			{
			cout << "Simulacion en t = " << (delta_t * s) << endl;
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
			// Calcula la fuerza total sobre q
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
			// Calcula la posicion y velocidad de q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Resultado final\n";
	printPosAndVel(N, pos, vel);

	for(int i=0; i<N; i++)
		omp_destroy_lock(&locks[i]);
}

void nBody_parallel_reduced_simple(string fileName, int N, int delta_t, int numSim, int tipoRes)
{
	// Creando el vector de masas
	double *masses;
	masses = new double[N];

	// Creando los vectores para almacenar las posiciones y las velocidades
	double **pos, **vel, **forces;
	pos = new double *[N];
	vel = new double *[N];
	forces = new double*[N];

	for(int i = 0; i < N; i++)
	{
		pos[i] = new double[2];
		vel[i] = new double[2];
		forces[i] = new double[2];
		forces[i][X] = forces[i][Y] = 0;
	}

	if(!loadDataFromFile(fileName, N, masses, pos, vel))
		return;

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

#pragma omp parallel
	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
#pragma omp single
			{
			cout << "Simulacion en t = " << (delta_t * s) << endl;
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
#pragma omp for
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
}

void nBody_parallel_reduced_critical(string fileName, int N, int delta_t, int numSim, int tipoRes)
{
	// Creando el vector de masas
	double *masses;
	masses = new double[N];

	// Creando los vectores para almacenar las posiciones y las velocidades
	double **pos, **vel, **forces;
	pos = new double *[N];
	vel = new double *[N];
	forces = new double*[N];

	for(int i = 0; i < N; i++)
	{
		pos[i] = new double[2];
		vel[i] = new double[2];
		forces[i] = new double[2];
		forces[i][X] = forces[i][Y] = 0;
	}

	if(!loadDataFromFile(fileName, N, masses, pos, vel))
		return;

	double x_diff, y_diff, dist, dist_cubed;
	double force_qkX, force_qkY;

#pragma omp parallel
	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
#pragma omp single
			{
			cout << "Simulacion en t = " << (delta_t * s) << endl;
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
			// Calcula la fuerza total sobre q
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
			// Calcula la posicion y velocidad de q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Resultado final\n";
	printPosAndVel(N, pos, vel);
}

void nBody_parallel_reduced_locforces(string fileName, int N, int delta_t, int numSim, int tipoRes)
{
	// Creando el vector de masas
	double *masses;
	masses = new double[N];

	// Creando los vectores para almacenar las posiciones y las velocidades
	double **pos, **vel, **forces;
	pos = new double *[N];
	vel = new double *[N];
	forces = new double*[N];

	for(int i = 0; i < N; i++)
	{
		pos[i] = new double[2];
		vel[i] = new double[2];
		forces[i] = new double[2];
		forces[i][X] = forces[i][Y] = 0;
	}

	if(!loadDataFromFile(fileName, N, masses, pos, vel))
		return;

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
				cout << "Simulacion en t = " << (delta_t * s) << endl;
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

			// Calcula la fuerza total sobre q
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
			// Calcula la posicion y velocidad de q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}
	}
	cout << "Resultado final\n";
	printPosAndVel(N, pos, vel);
}

int main() {
	int N, delta_t, numSim, tipoRes;
	string fileName = "";

	cout << "Ingrese los siguientes datos:";
	cout << "\n- Numero de particulas: ";
	cin >> N;
	cout << "\n- Intervalo de tiempo (s): ";
	cin >> delta_t;
	cout << "\n- Numero de iteraciones: ";
	cin >> numSim;
	cout << "\n- Archivo con los datos: ";
	cin >> fileName;
	cout << "\nResultado final (0) / Resultado por iteracion (1): ";
	cin >> tipoRes;

	//nBody_parallel_reduced_simple(fileName, N, delta_t, numSim, tipoRes);
	//nBody_parallel_reduced_critical(fileName, N, delta_t, numSim, tipoRes);
	//nBody_parallel_reduced_locks(fileName, N, delta_t, numSim, tipoRes);
	nBody_parallel_reduced_locforces(fileName, N, delta_t, numSim, tipoRes);
}
