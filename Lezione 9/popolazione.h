#ifndef __popolazione__
#define __popolazione__

#include <armadillo>
#include "random.h"
#include "percorso.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using namespace arma;

// Class that manages a population of TSP routes for a Genetic Algorithm
class Population {

private:
    Random _rnd;                      // Random number generator
    int _npop;                        // Number of individuals in the population
    const double _selexp = 2;         // Selection pressure exponent
    arma::field<Route> _pop;          // Container for the population of routes
    arma::mat _dist_matrix;           // Distance matrix between cities
    const double _pcross = 0.90;      // Crossover probability
    int _ngen;                        // Number of generations to evolve

public:
    void initialize_pop(int npop, int ngen, arma::mat* dist_matrix); // Initialize population
    Route get_percorso(int i);                                       // Get i-th route
    int size() const;                                                // Return population size
    void sort_by_length();                                           // Sort population by route length
    int select_index();                                              // Return index based on selection probability
    Route select();                                                  // Select a route using selection operator
    arma::field<Route> crossover(Route a, Route b);                  // Perform crossover between two parents
    arma::Col<int> sort_by_reference(arma::Col<int> a, arma::Col<int> ref); // Order vector 'a' based on order in 'ref'
    void evolve(const arma::mat dist_matrix);                        // Run the Genetic Algorithm
    int pbc(int city) const;                                         // Periodic boundary condition for index wrapping
};

#endif
