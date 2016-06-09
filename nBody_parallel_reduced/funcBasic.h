#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

#define DIM 2
#define X 	0
#define Y 	1

void initData(int num, double *masses, double pos[][2], double vel[][2])
{
	for(int i=0; i<num; i++)
	{ // (max - min) * ( (double)rand() / (double)RAND_MAX ) + min
		masses[i] = (1000 - 900) * rand() + 900;
		pos[i][X] = (400 - 200) * rand() + 200;
		pos[i][Y] = (400 - 200) * rand() + 200;
		vel[i][X] = (500 - 200) * rand() + 200;
		vel[i][Y] = (500 - 200) * rand() + 200;
	}
}

bool loadDataFromFile(string fileName, int num, double *mass, double pos[][2], double vel[][2])
{
	fstream myFile(fileName.c_str(), std::ios_base::in);
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

void printPosAndVel(int num, double pos[][2], double vel[][2])
{
	cout << "Particle \t Position \t Velocity\n";
	for(int i=0; i<num; i++)
	{
		cout << (i+1) << "\t(" << pos[i][X] << "," << pos[i][Y] << ")\t "
			<< "(" << vel[i][X] << "," << vel[i][Y] << ")\n";
	}
}

void clearForces(int num, double forc[][2])
{
	for(int i = 0; i < num; i++)
	{
		forc[i][X] = forc[i][Y] = 0;
	}
}
