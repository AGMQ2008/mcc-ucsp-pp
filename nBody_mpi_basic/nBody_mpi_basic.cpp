#include <iostream>
#include <mpi.h>
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
	int my_rank, comm_sz;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	MPI_Datatype vect_mpi_t;
	MPI_Type_contiguous(DIM, MPI_DOUBLE, &vect_mpi_t);
	MPI_Type_commit(&vect_mpi_t);

	MPI_Comm comm = MPI_COMM_WORLD;
	int loc_n = N / comm_sz;

	double loc_vel[1][2];
	loc_vel[1][X] = 0;
	loc_vel[1][Y] = 0;

	/*
	loc_vel = new double *[loc_n];

	for(int i=0; i<loc_n; i++)
	{
		loc_vel[i] = new double[2];
		loc_vel[i][X] = 0;
		loc_vel[i][Y] = 0;
	}
	*/

	if(my_rank == 0)
	{cout << "master executed";
		// Loading initial data
		initData(N, masses, pos, vel);
		// Clear forces array
		clearForces(N, forces);
	}

	int loc_pos = my_rank * loc_n;

	MPI_Bcast(masses, N, MPI_DOUBLE, 0, comm);
	MPI_Bcast(pos, N, vect_mpi_t, 0, comm);
	MPI_Bcast(vel, N, vect_mpi_t, 0, comm);
	MPI_Scatter(vel, loc_n, vect_mpi_t, loc_vel, loc_n, vect_mpi_t, 0, comm);

	double x_diff, y_diff, dist, dist_cubed;

	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
			MPI_Gather(loc_vel, loc_n, vect_mpi_t, vel, loc_n, vect_mpi_t, 0, comm);
			//MPI_Allgather(MPI_IN_PLACE, loc_n, vect_mpi_t, vel, loc_n, vect_mpi_t, comm);
			if (my_rank == 0)
			{
				cout << "Simulacion en t = " << (delta_t * s) << endl;
				printPosAndVel(N, pos, vel);
			}
		}

		for(int q=loc_pos; q<(loc_pos + loc_n); q++)
		{
			forces[q][X] = forces[q][Y] = 0;

			// Calcula la fuerza total sobre q
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
		int i=0;
		for(int q=loc_pos; q<(loc_pos+loc_n); q++)
		{
			// Calcula la posicion y velocidad de q
			pos[q][X] += delta_t * loc_vel[i][X];
			pos[q][Y] += delta_t * loc_vel[i][Y];
			loc_vel[i][X] += delta_t / masses[q] * forces[q][X];
			loc_vel[i][Y] += delta_t / masses[q] * forces[q][Y];
			i++;
		}

		MPI_Allgather(MPI_IN_PLACE, loc_n, vect_mpi_t, pos, loc_n, vect_mpi_t, comm);
		//MPI_Allgather(MPI_IN_PLACE, loc_n, vect_mpi_t, vel, loc_n, vect_mpi_t, comm);
	}

	MPI_Gather(loc_vel, loc_n, vect_mpi_t, vel, loc_n, vect_mpi_t, 0, comm);
	if(my_rank == 0)
	{
		cout << "Resultado final\n";
		printPosAndVel(N, pos, vel);
	}

	MPI_Finalize();

	return 0;
}
