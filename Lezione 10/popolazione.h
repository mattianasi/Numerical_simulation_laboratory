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
    Random _rnd;
    int _npop;
    const double _selexp = 2;
    arma::field <Route> _pop;
    arma::mat _dist_matrix;
    const double _pcross = 0.90;
    int _ngen;

    public:
    void initialize_pop(int npop, int ngen, arma::mat * dist_matrix);
    Route get_percorso(int i);
    int size() const;
    void sort_by_length();
    int select_index();
    Route select();
    arma::field<Route> crossover(Route a, Route b);
    arma::Col<int> sort_by_reference(arma::Col<int> a, arma::Col<int> ref);
    void evolve(const arma::mat dist_matrix, int rank, int size, MPI_Comm comm, int Nmigr) ;
    int pbc(int city) const;
    void replace_best(const Route& migrant);
    void migrate(int rank, int size, MPI_Comm comm);
    double get_best_length() const {
        return _pop[0].calculate_length(&_dist_matrix);
    }
    arma::Col<int> get_best_route()  {
        return _pop[0].getroute();
    }

};

#endif