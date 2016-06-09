#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

#define DIM 2
#define X 	0
#define Y 	1

double genNumRandom(int v_min, int v_max)
{
	return (v_max - v_min) * ((double)rand() / (double)RAND_MAX) + v_min;
}

void initData(int num, double *masses, double pos[][DIM], double vel[][DIM])
{
	for(int i=0; i<num; i++)
	{
		masses[i] = genNumRandom(900, 1000);
		pos[i][X] = genNumRandom(200, 400);
		pos[i][Y] = genNumRandom(200, 400);
		vel[i][X] = genNumRandom(200, 500);
		vel[i][Y] = genNumRandom(200, 500);
	}
}

bool loadDataFromFile(string fileName, int num, double *mass, double pos[][DIM], double vel[][DIM])
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

void copyVector(int sz, double source[][DIM], double dest[][DIM])
{
	for(int i = 0; i < sz; i++)
	{
		dest[i][X] = source[i][X];
		dest[i][Y] = source[i][Y];
	}
}

void clearVector(int sz, double v[][DIM])
{
	for(int i = 0; i < sz; i++)
	{
		v[i][X] = v[i][Y] = 0;
	}
}
