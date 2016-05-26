#include <iostream>
#include <mpi.h>
#include <fstream>
#include <string>

using namespace std;

#define DIM 2
#define X 0
#define Y 1

double const G = 6.673;


void printPosAndVel(int N, double** pos, double** vel)
{
	cout << "Particle \t Position \t Velocity\n";
	for(int i=0; i<N; i++)
	{
		cout << (i+1) << "\t(" << pos[i][X] << "," << pos[i][Y] << ")\t "
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

int main(int argc, char *argv[]) {
	int my_rank, size;

	double *masses, **pos, **vel, **forces;

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

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(my_rank == 0)
	{
		masses = new double[N];
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
			return 0;
	}

	MPI_Comm comm = MPI_COMM_WORLD;
	int loc_n = N / size;
	int loc_vel = N / size;
	MPI_Datatype vect_mpi_t;
	MPI_Type_contiguous(DIM, MPI_DOUBLE, &vect_mpi_t);
	MPI_Type_commit(&vect_mpi_t);

	MPI_Bcast(masses, N, MPI_DOUBLE, 0, comm);
	MPI_Bcast(pos, N, vect_mpi_t, 0, comm);
	MPI_Scatter(vel, loc_n, vect_mpi_t, vel, loc_n, vect_mpi_t, 0, comm);

	double x_diff, y_diff, dist, dist_cubed;

	for(int s=0; s<numSim; s++)
	{
		if (tipoRes == 1)
		{
			if (my_rank == 0)
			{
				cout << "Simulacion en t = " << (delta_t * s) << endl;
				printPosAndVel(N, pos, vel);
			}
		}

		for(int q=0; q<N; q++)
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
		for(int q=0; q<N; q++)
		{
			// Calcula la posicion y velocidad de q
			pos[q][X] += delta_t * vel[q][X];
			pos[q][Y] += delta_t * vel[q][Y];
			vel[q][X] += delta_t / masses[q] * forces[q][X];
			vel[q][Y] += delta_t / masses[q] * forces[q][Y];
		}

		MPI_Allgather(MPI_IN_PLACE, loc_n, vect_mpi_t, pos, loc_n, vect_mpi_t, comm);
	}

	if(my_rank == 0)
	{
		cout << "Resultado final\n";
		printPosAndVel(N, pos, vel);
	}

	MPI_Finalize();

	return 0;
}
