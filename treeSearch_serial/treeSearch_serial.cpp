#include <iostream>
#include <stdlib.h>

using namespace std;


// Global variables
const int N = 4;	// Number of cities
const int hometown = 0;	// vertex o city 0, salesperson's hometown
int **digraph;	// Data structure of input digraph

typedef int city_t;

typedef struct
{
	city_t cities[N];
	int nCities;
	double vPath;
} tour_t;

tour_t best_tour;	// best tour with minimal cost

void createDigraph()
{
	digraph = new int*[N];
	for(int i=0; i<N; i++)
	{
		digraph[i] = new int[N];
		for(int j=0; j<N; j++)
		{
			if(i == j)
				digraph[i][j] = 0;
			else
				digraph[i][j] = (rand() % 5) + 1;
		}
	}
}

// Examines if there are n cities on the partial tour
int City_count(tour_t tour)
{
	return tour.nCities;
}

bool Best_tour(tour_t tour)
{
	return tour.vPath < best_tour.vPath;
}

void Update_best_tour(tour_t tour)
{
	best_tour.nCities = tour.nCities;
	best_tour.vPath = tour.vPath;
	for(int i=0; i<tour.nCities; i++)
		best_tour.cities[i] = tour.cities[i];
}

// Checks if the city has already been visited
bool Feasible(tour_t tour, city_t city)
{
	for(int i=0; i<tour.nCities; i++)
	{
		if(tour.cities[i] == city)
			return false;
	}
	return true;
}

void Add_city(tour_t *tour, city_t city)
{
	int lastCity = tour->cities[tour->nCities-1];
	tour->cities[tour->nCities++] = city;
	tour->vPath += digraph[lastCity][city];
}

void Remove_last_city(tour_t *tour, city_t city)
{
	int lastCity = tour->cities[tour->nCities-2];
	tour->cities[tour->nCities--] = -1;
	tour->vPath -= digraph[lastCity][city];
}

void Depth_first_search(tour_t tour)
{
	city_t city;

	if(City_count(tour) == N)
	{
		if(Best_tour(tour))
			Update_best_tour(tour);
	}
	else
	{
		// Se visita a los vecinos
		for(city=0; city<N; city++)
		{
			if(Feasible(tour, city))
			{
				Add_city(&tour, city);
				Depth_first_search(tour);
				Remove_last_city(&tour, city);
			}
		}
	}
}

int main()
{
	createDigraph();

	cout << "El grafo generado de " << N << " ciudades es:"<<endl;
	for(int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
			cout << digraph[i][j] << " ";
		cout << endl;
	}

	//best_tour.cities = new city_t[N];
	best_tour.nCities = 0;
	best_tour.vPath = 9999999;

	cout << "seteado best tour";
	tour_t tour;
	//tour.cities = new city_t[N];
	tour.nCities = 0;
	tour.vPath = 0;

	cout << "Se agregara la ciudad de inicio";
	Add_city(&tour, hometown);

	Depth_first_search(tour);

	for(int i=0; i<best_tour.nCities; i++)
	{
		cout << best_tour.cities[i] << " - ";
	}
	cout << "Longitud: " << best_tour.vPath;

	return 0;
}
