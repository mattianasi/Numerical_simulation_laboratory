#include <armadillo>
#include "random.h"
#include "percorso.h"
#include "popolazione.h"
#include <iostream>
#include <cmath>
#include <algorithm> 


// Initialize the population with npop individuals, to evolve for ngen generations
void Population :: initialize_pop(int npop, int ngen, arma::mat * dist_matrix){
    _npop=npop, _ngen=ngen;

    int p1, p2; // Prime numbers used to initialize the RNG
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    int seed[4]; // Seed for the RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed,p1,p2);

    _dist_matrix = *dist_matrix; // Store the distance matrix

    _pop.set_size(_npop);
    for (int i = 0; i < _npop; ++i) {
        _pop[i].initialize(&_dist_matrix, _rnd); // Create and initialize each individual
    }
}

// Return the i-th individual from the population
Route Population :: get_percorso(int i){
    return _pop(i);
}

// Select an individual index based on power-law probability
int Population::select_index() {
    double r = _rnd.Rannyu(); 
    int j = int(_npop * std::pow(r, _selexp)); // Biases toward lower indices (better individuals)
    return j;
}

// Sort population by tour length in ascending order
void Population::sort_by_length() {
    std::vector<Route> temp_vec(_npop);
    for (int i = 0; i < _npop; ++i){
        temp_vec[i] = _pop[i];
    }

    // Compare based on calculated route length
    std::sort(temp_vec.begin(), temp_vec.end(), [&](const Route& a, const Route& b) {
        return a.calculate_length(&_dist_matrix) < b.calculate_length(&_dist_matrix);
    });
    
    for (int i = 0; i < _npop; ++i){
        _pop[i] = temp_vec[i];
    }
}

// Select and return one individual based on selection strategy
Route Population::select() { 
    int j = select_index();
    return _pop(j);
}

// Return current population size
int Population::size() const {
    return _pop.size();
}

// Periodic boundary condition for cities index wrapping
int Population :: pbc(int city) const{
    int ndim = 34;
    if(city >= ndim){
        return city - ndim + 1;
    }
    return city;
}

// Partially matched crossover between two parent routes
arma::field<Route> Population::crossover(Route a, Route b) {
    if (_rnd.Rannyu() < _pcross) {
        int n_cities = a.getdim();
        arma::Col<int> parent1 = a.getroute();
        arma::Col<int> parent2 = b.getroute();

        int cut1 = int(_rnd.Rannyu(1, n_cities - 2)); // Ensure cuts are internal
        int cut2 = int(_rnd.Rannyu(cut1 + 1, n_cities)); // Second cut after the first

        arma::Col<int> offspring1(n_cities, arma::fill::zeros);
        arma::Col<int> offspring2(n_cities, arma::fill::zeros);

        // Copy segment between cuts
        for (int i = cut1; i <= cut2; ++i) {
            offspring1[i] = parent1[i];
            offspring2[i] = parent2[i];
        }

        // Fill remaining cities while avoiding duplicates
        int idx1 = (cut2 + 1) % n_cities;
        int idx2 = idx1;
        for (int i = 0; i < n_cities; ++i) {
            int city = parent2[(cut2 + 1 + i) % n_cities];
            if (arma::any(offspring1 == city)) continue;
            offspring1[idx1] = city;
            idx1 = (idx1 + 1) % n_cities;
        }

        for (int i = 0; i < n_cities; ++i) {
            int city = parent1[(cut2 + 1 + i) % n_cities];
            if (arma::any(offspring2 == city)) continue;
            offspring2[idx2] = city;
            idx2 = (idx2 + 1) % n_cities;
        }

        arma::field<Route> children(2);
        children[0].initialize(&_dist_matrix, _rnd);
        children[1].initialize(&_dist_matrix, _rnd);
        children[0].setroute(offspring1);
        children[1].setroute(offspring2);

        return children;
    } else {
        return {a, b}; // Skip crossover, return unmodified parents
    }
}

// Reorders vector `a` based on the relative order of `ref`
arma::Col<int> Population :: sort_by_reference(arma::Col<int> a, arma::Col<int> ref){
    arma::Col<int> ref_pos(ref.size()); 
    for(int i = 0; i < ref.size(); i++){
        ref_pos[ref[i]-1] = i;
    }

    vector<int> sorted(a.begin(),a.end());
    sort(sorted.begin(), sorted.end(), [&ref_pos](int j, int k){
        return ref_pos[j] < ref_pos[k];
    });

    return arma::Col<int>(sorted); // Convert back to Armadillo column
}

// Main evolution loop of the genetic algorithm
void Population::evolve(const arma::mat dist_matrix) {
    ofstream coutf("results.dat");
    coutf << "#" << "\t\t L1 \t\t <L1>" << endl;

    sort_by_length(); // Initial sort

    for (int i = 0; i < _ngen; i++){
        arma::field<Route> figli(_npop);

        // Generate new individuals via selection and crossover
        for(int j = 0; j < _npop/2; j++){
            Route a = select();
            Route b = select();
            arma::field<Route> temp = crossover(a,b);
            figli[j].initialize(&_dist_matrix, _rnd);
            figli[j+_npop/2].initialize(&_dist_matrix, _rnd);
            figli[j].setroute(temp[0].getroute());
            figli[j+_npop/2].setroute(temp[1].getroute());
        }

        // Apply mutation to each offspring
        for(int j = 0; j < _npop; j++){
            figli[j].mutate();
        }

        _pop = figli;
        sort_by_length(); // Sort population after mutation

        // Track statistics
        double sum = 0;
        for (int j = 0; j < _npop/2; j++){
            sum += _pop[j].calculate_length(&dist_matrix);
        }
        coutf << i << "\t\t" << _pop[0].calculate_length(&dist_matrix) << "\t\t" << sum/(_npop/2) << endl;
    }

    coutf.close();

    // Output best route to file
    ofstream out("best_route.dat");
    out << "#Best route \t \t x \t \t y" << endl;
    for (int i = 0; i < _pop[0].getdim(); i++){
        out << _pop[0].getroute()[i] << endl;
    }
    out.close();
}
