#include <iostream>
#include <stdlib.h>
#include <stack>
#include <time.h>

using namespace std;

#define N 10	// Number of cities
#define hometown 0		// salesperson's hometown (city 0)
#define MAX_DIST 9999999

typedef struct
{
	int cities[N];
	int nCities;
	double vPath;
} tour_t;

// Global variables
int digraph[N][N];		// Data structure of input digraph
tour_t best_tour;	// best tour with minimal cost

#define Tour_city(tour, i) (tour->cities[i])


// This function creates a digraph with random values
void createDigraph()
{
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			if(i == j)
				digraph[i][j] = 0;
			else
				digraph[i][j] = (rand() % 5) + 1;
		}
	}
}

// This function prints the cities and length of the path
void printPath(tour_t tour)
{
	cout << "Path: ";
	for(int i=0; i<tour.nCities; i++)
	{
		cout << (tour.cities[i]) << " - ";
	}
	cout << tour.cities[0];
	cout << "\nLength: " << tour.vPath << endl;
}

// This function prints the values of digraph in console
void printDigraph()
{
	cout << "Digraph with " << N << " cities:" << endl;
	for(int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
			cout << digraph[i][j] << " ";
		cout << endl;
	}
}

// Return the number of cities in the path or tour
int City_count(tour_t tour)
{
	return tour.nCities;
}

// This function verifies if the tour is better than the current best tour
bool Is_best_tour(tour_t tour)
{
	return tour.vPath < best_tour.vPath;
}

// Updates the current best tour to tour
void Update_best_tour(tour_t tour)
{
	best_tour.nCities = tour.nCities;
	best_tour.vPath = tour.vPath;
	for(int i=0; i<tour.nCities; i++)
		best_tour.cities[i] = tour.cities[i];
}

// Checks if the city has already been visited
bool Feasible(tour_t tour, int city)
{
	for(int i=0; i<tour.nCities; i++)
	{
		if(tour.cities[i] == city)
			return false;
	}
	return true;
}

// Adds a new city to the tour
void Add_city(tour_t *tour, int city)
{
	int iniCity, lastCity;

	if(tour->nCities > 0)
	{
		iniCity = Tour_city(tour, 0);
		lastCity = Tour_city(tour, tour->nCities - 1);
		tour->vPath = tour->vPath + digraph[lastCity][city] + digraph[city][iniCity] - digraph[lastCity][iniCity];
	}
	tour->cities[tour->nCities++] = city;
}

// Removes the last city in the tour
void Remove_last_city(tour_t *tour)
{
	if(tour->nCities <= 1)
	{
		tour->nCities = 0;
		tour->vPath = 0;
		return;
	}

	int iniCity, lastCity, newLastCity;

	iniCity = Tour_city(tour, 0);
	lastCity = Tour_city(tour, tour->nCities - 1);
	newLastCity = Tour_city(tour, tour->nCities - 2);

	tour->nCities--;
	tour->vPath = tour->vPath - (digraph[newLastCity][lastCity] + digraph[lastCity][iniCity]) + digraph[newLastCity][iniCity];
}

// Recursive DFS
void Depth_first_search(tour_t tour)
{
	if(City_count(tour) == N)
	{
		if(Is_best_tour(tour))
			Update_best_tour(tour);
	}
	else
	{
		// Visiting neighbors
		for(int city=0; city<N; city++)
		{
			if(Feasible(tour, city))
			{
				Add_city(&tour, city);
				Depth_first_search(tour);
				Remove_last_city(&tour);
			}
		}
	}
}

#define NO_CITY -1

// Iterative DFS using the variable NO_CITY
void Depth_first_search_iter01(tour_t tour)
{
	stack<int> myStack;

	for(int city = N-1; city >= 1; city--)
	{
		myStack.push(city);

		while(!myStack.empty())
		{
			city = myStack.top();
			myStack.pop();

			if(city == NO_CITY)
			{
				Remove_last_city(&tour);
			}
			else
			{
				Add_city(&tour, city);
				if(City_count(tour) == N)
				{
					if(Is_best_tour(tour))
						Update_best_tour(tour);
					Remove_last_city(&tour);
				}
				else
				{
					myStack.push(NO_CITY);
					for(int nbr = N-1; nbr>=1; nbr--)
					{
						if(Feasible(tour, nbr))
							myStack.push(nbr);
					}
				}
			}
		}
	}
}

// Iterative DFS without NO_CITY
void Depth_first_search_iter02(tour_t tour)
{
	stack<tour_t> myStack;
	tour_t curr_tour;

	myStack.push(tour);

	while(!myStack.empty())
	{
		curr_tour = myStack.top();
		myStack.pop();

		if(City_count(curr_tour) == N)
		{
			if(Is_best_tour(curr_tour))
				Update_best_tour(curr_tour);
		}
		else
		{
			for(int nbr = N-1; nbr >= 1; nbr--)
			{
				if(Feasible(curr_tour, nbr))
				{
					Add_city(&curr_tour, nbr);
					myStack.push(curr_tour);
					Remove_last_city(&curr_tour);
				}
			}
		}
	}
}

int main()
{
	best_tour.nCities = 0;
	best_tour.vPath = MAX_DIST;

	createDigraph();
	//printDigraph();

	tour_t tour;
	tour.nCities = 0;
	tour.vPath = 0;
	Add_city(&tour, hometown);

	clock_t begin, end;
	double elapsed_secs = 0;
	begin = clock();

	Depth_first_search(tour);
	//Depth_first_search_iter01(tour);
	//Depth_first_search_iter02(tour);

	end = clock();
	elapsed_secs = double(end - begin) * 1.0 / CLOCKS_PER_SEC;

	printPath(best_tour);
	cout << "Time elapsed: " << elapsed_secs << " secs";

	return 0;
}
