#include <iostream>
#include <stdlib.h>
#include <stack>
#include <queue>
#include <pthread.h>

using namespace std;

#define N 6		// Number of cities
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
pthread_mutex_t best_tour_lock = PTHREAD_MUTEX_INITIALIZER;

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
void *Depth_first_search_static(void *t)
{
	tour_t *tour;
	tour = (tour_t *) t;

	stack<tour_t> my_stack;
	tour_t curr_tour;
	my_stack.push(*tour);

	while(!my_stack.empty())
	{
		curr_tour = my_stack.top();
		my_stack.pop();

		if(City_count(curr_tour) == N)
		{
			pthread_mutex_lock(&best_tour_lock);
			if(Is_best_tour(curr_tour))
				Update_best_tour(curr_tour);
			pthread_mutex_unlock(&best_tour_lock);
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

int threads_in_cond_wait = 0;
pthread_mutex_t term_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t term_cond_var = PTHREAD_COND_INITIALIZER;
stack<tour_t> new_stack;

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
	if(my_stack->size() >= 2 && threads_in_cond_wait > 0 && new_stack.size() == 0)
	{
		pthread_mutex_lock(&term_mutex);
		if(threads_in_cond_wait > 0 && new_stack.size() == 0)
		{
			// split my_stack creating new_stack
			int sStack = my_stack->size();
			for(int i=0; i<sStack/2; i++)
			{
				new_stack.push(my_stack->top());
				my_stack->pop();
			}
			pthread_cond_signal(&term_cond_var);
		}
		pthread_mutex_unlock(&term_mutex);
		return 0;
	}
	else if(!my_stack->empty())
		return 0;
	else // stack is empty
	{
		pthread_mutex_lock(&term_mutex);
		if(threads_in_cond_wait == thread_count - 1)
		{
			threads_in_cond_wait++;
			pthread_cond_broadcast(&term_cond_var);
			pthread_mutex_unlock(&term_mutex);
			return 1;
		}
		else
		{
			threads_in_cond_wait++;
			while(pthread_cond_wait(&term_cond_var, &term_mutex) != 0)
				;

			if(threads_in_cond_wait < thread_count)
			{
				setStack(my_stack, new_stack);
				clearNewStack();
				threads_in_cond_wait--;
				pthread_mutex_unlock(&term_mutex);
				return 0;
			}
			else
			{
				pthread_mutex_unlock(&term_mutex);
				return 1;
			}
		}
	}
}

void *Depth_first_search_dynamic(void *t)
{
	tour_t *tour;
	tour = (tour_t *) t;

	stack<tour_t> my_stack;
	tour_t curr_tour;
	my_stack.push(*tour);

	while(!Terminated(&my_stack))
	{
		curr_tour = my_stack.top();
		my_stack.pop();

		if(City_count(curr_tour) == N)
		{
			pthread_mutex_lock(&best_tour_lock);
			if(Is_best_tour(curr_tour))
				Update_best_tour(curr_tour);
			pthread_mutex_unlock(&best_tour_lock);
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
	best_tour.nCities = 0;
	best_tour.vPath = MAX_DIST;

	createDigraph();
	//printDigraph();

	tour_t tour;
	tour.nCities = 0;
	tour.vPath = 0;
	Add_city(&tour, hometown);

	thread_count = BFS(tour);

	// Creation of the threads, one thread per tour in the initial queue
	pthread_t threads[thread_count];
	for(int i=0; i<thread_count; i++)
	{
		tour = myQueue.front();
		//pthread_create(&threads[i], NULL, Depth_first_search_static, (void *)&tour);
		pthread_create(&threads[i], NULL, Depth_first_search_dynamic, (void *)&tour);
		myQueue.pop();
	}

	for(int i=0; i<thread_count; i++)
		pthread_join(threads[i], NULL);

	printPath(best_tour);

	pthread_mutex_destroy(&best_tour_lock);
	pthread_mutex_destroy(&term_mutex);
	pthread_cond_destroy(&term_cond_var);

	return 0;
}
