#include <iostream>
#include <stdlib.h>
#include <stack>
#include <queue>
#include <omp.h>

using namespace std;

#define N 10		// Number of cities
#define hometown 0		// salesperson's hometown (city 0)
#define MAX_DIST 9999999
#define NUM_MIN_THREADS 6

typedef struct
{
	int cities[N];
	int nCities;
	double vPath;
} tour_t;

// Global variables
int digraph[N][N];	// Data structure of input digraph
tour_t best_tour;	// best tour with minimal cost

int thread_count = 0;
queue<tour_t> myQueue;
omp_lock_t best_tour_lock;

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

// breadth-first-search returns the number of tours
int BFS(tour_t tour)
{
	myQueue.push(tour);

	while(!myQueue.empty())
	{
		tour = myQueue.front();
		myQueue.pop();
		if(City_count(tour) == N)
		{
			if(Is_best_tour(tour))
				Update_best_tour(tour);
		}
		else
		{
			for(int city = 0; city < N; city++)
			{
				if(Feasible(tour, city))
				{
					Add_city(&tour, city);
					myQueue.push(tour);
					Remove_last_city(&tour);
				}
			}
			if(myQueue.size() >= NUM_MIN_THREADS)
				break;
		}
	}
	return myQueue.size();
}

// Iterative function used by the threads to find the best tour
void Depth_first_search_static(tour_t tour)
{
	stack<tour_t> my_stack;
	tour_t curr_tour;
	my_stack.push(tour);

	while(!my_stack.empty())
	{
		curr_tour = my_stack.top();
		my_stack.pop();

		if(City_count(curr_tour) == N)
		{
#pragma omp critical
			{
				if(Is_best_tour(curr_tour))
					Update_best_tour(curr_tour);
			}
		}
		else
		{
			for(int city = N-1; city >= 1; city--)
			{
				if(Feasible(curr_tour, city))
				{
					Add_city(&curr_tour, city);
					my_stack.push(curr_tour);
					Remove_last_city(&curr_tour);
				}
			}
		}
	}
}

stack<tour_t> new_stack;
queue<int> threadsSleep;

int awakened_thread = -1;
int work_remains = 1;
omp_lock_t term_lock;

void clearNewStack()
{
	int size = new_stack.size();
	for(int i=0; i<size; i++)
		new_stack.pop();
}

void setStack(stack<tour_t> *s, stack<tour_t> d)
{
	stack<tour_t> tmp;
	int size = d.size();
	for(int i=0; i<size; i++)
		tmp.push(d.top());

	for(int i=0; i<size; i++)
	{
		s->push(tmp.top());
	}
}

bool Terminated(stack<tour_t> *my_stack)
{
	int got_lock;
	int my_rank = omp_get_thread_num();

	if(my_stack->size() >= 2 && threadsSleep.size() > 0 && new_stack.size() == 0)
	{
		got_lock = omp_test_lock(&term_lock);
		if(got_lock != 0)
		{
			if(threadsSleep.size() > 0 && new_stack.size() == 0)
			{
				// split my_stack creating new_stack
				int sStack = my_stack->size();
				for(int i=0; i<sStack/2; i++)
				{
					new_stack.push(my_stack->top());
					my_stack->pop();
				}
				awakened_thread = threadsSleep.front();
			}
			omp_unset_lock(&term_lock);
			return 0;
		}
	}
	else if(!my_stack->empty())
		return 0;
	else // stack is empty
	{
		omp_set_lock(&term_lock);
		if(threadsSleep.size() == (thread_count - 1))
		{
			threadsSleep.push(my_rank);
			work_remains = 0;
			omp_unset_lock(&term_lock);
			return 1;
		}
		else
		{
			threadsSleep.push(my_rank);

			omp_unset_lock(&term_lock);
			while(awakened_thread != my_rank && work_remains)
				;
			omp_set_lock(&term_lock);

			if(threadsSleep.size() < thread_count)
			{
				setStack(my_stack, new_stack);
				clearNewStack();
				threadsSleep.pop();
				awakened_thread = -1;
				omp_unset_lock(&term_lock);
				return 0;
			}
			else
			{
				omp_unset_lock(&term_lock);
				return 1;
			}
		}
	}
}

void Depth_first_search_dynamic(tour_t tour)
{
	stack<tour_t> my_stack;
	tour_t curr_tour;
	my_stack.push(tour);

	while(!Terminated(&my_stack))
	{
		curr_tour = my_stack.top();
		my_stack.pop();

		if(City_count(curr_tour) == N)
		{
			omp_set_lock(&best_tour_lock);
			if(Is_best_tour(curr_tour))
				Update_best_tour(curr_tour);
			omp_unset_lock(&best_tour_lock);
		}
		else
		{
			for(int city = N-1; city >= 1; city--)
			{
				if(Feasible(curr_tour, city))
				{
					Add_city(&curr_tour, city);
					my_stack.push(curr_tour);
					Remove_last_city(&curr_tour);
				}
			}
		}
	}
}

int main()
{
	tour_t tour;

	best_tour.nCities = 0;
	best_tour.vPath = MAX_DIST;

	createDigraph();
	//printDigraph();

	tour.nCities = 0;
	tour.vPath = 0;
	Add_city(&tour, hometown);

	thread_count = BFS(tour);
	cout << "Number of threads: " << thread_count << endl;

	omp_set_num_threads(thread_count);
	omp_init_lock(&best_tour_lock);

	double start, end;
	start = omp_get_wtime();

#pragma omp parallel for
	for(int i = 0; i < thread_count; i++)
	{
		tour = myQueue.front();
		//Depth_first_search_static(tour);
		Depth_first_search_dynamic(tour);
		myQueue.pop();
	}

	end = omp_get_wtime();

#pragma omp single
	{
		printPath(best_tour);
		cout << "Time elapsed: " << (end - start) << endl;
	}

	return 0;
}
