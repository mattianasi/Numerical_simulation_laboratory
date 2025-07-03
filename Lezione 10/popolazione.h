#ifndef __popolazione__
#define __popolazione__

#include <armadillo>
#include "random.h"
#include "percorso.h"
#include <iostream>
#include <vector>
#include <mpi.h>
#include <cmath>

using namespace std;
using namespace arma;

class Population {

    private:
    Random _rnd;                            // Random number generator used by the population
    int _npop;                              // Number of individuals in the population
    const double _selexp = 2;               // Selection exponent for biased selection
    arma::field <Route> _pop;              // Container holding all individuals (routes)
    arma::mat _dist_matrix;                // Distance matrix between cities
    const double _pcross = 0.90;           // Probability of performing crossover
    int _ngen;                              // Number of generations

    public:
    void initialize_pop(int npop, int ngen, arma::mat * dist_matrix); // Create initial population
    Route get_percorso(int i);             // Get individual i from the population
    int size() const;                      // Return population size
    void sort_by_length();                 // Sort individuals by tour length (fitness)
    int select_index();                    // Select index based on selection pressure
    Route select();                        // Return selected individual
    arma::field<Route> crossover(Route a, Route b); // Perform crossover between two parents
    arma::Col<int> sort_by_reference(arma::Col<int> a, arma::Col<int> ref); // Sort a by matching order in ref
    void evolve(const arma::mat dist_matrix, int rank, int size, MPI_Comm comm, int Nmigr); // Main GA evolution loop with MPI migration
    int pbc(int city) const;               // Apply periodic boundary conditions
    void replace_best(const Route& migrant); // Replace worst individual with a migrant (used in migration)
    void migrate(int rank, int size, MPI_Comm comm); // Perform migration between MPI processes
    double get_best_length() const {
        return _pop[0].calculate_length(&_dist_matrix); // Return best path length (after sorting)
    }
    arma::Col<int> get_best_route()  {
        return _pop[0].getroute(); // Return best route (after sorting)
    }

};

#endif
