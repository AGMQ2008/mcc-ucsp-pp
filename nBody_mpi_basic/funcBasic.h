#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

#define DIM 2
#define X 	0
#define Y 	1

void initData(int num, double *masses, double pos[][DIM], double vel[][DIM])
{
	for(int i=0; i<num; i++)
	{
		masses[i] = rand() % 3 + 1;
		pos[i][X] = rand() % 4;
		pos[i][Y] = rand() % 4;
		vel[i][X] = rand() % 5 + 1;
		vel[i][Y] = rand() % 5 + 1;
	}
}

bool loadDataFromFile(string fileName, int num, double *mass, double pos[][DIM], double vel[][DIM])
{
	fstream myFile(fileName, std::ios_base::in);
	double temp = 0.0;

	if(myFile.is_open())
	{
		for(int i = 0; i < num; i++)
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

void printPosAndVel(int num, double pos[][DIM], double vel[][DIM])
{
	cout << "Particle \t Position \t Velocity\n";
	for(int i=0; i<num; i++)
	{
		cout << (i+1) << "\t(" << pos[i][X] << "," << pos[i][Y] << ")\t "
			<< "(" << vel[i][X] << "," << vel[i][Y] << ")\n";
	}
}

void clearForces(int num, double forc[][DIM])
{
	for(int i = 0; i < num; i++)
	{
		forc[i][X] = forc[i][Y] = 0;
	}
}
