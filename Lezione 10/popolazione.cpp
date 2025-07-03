#include <armadillo>
#include "random.h"
#include "percorso.h"
#include "popolazione.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <mpi.h>
#include <algorithm> 


// Initialize the population of routes, RNG and distance matrix
void Population :: initialize_pop(int npop, int ngen, arma::mat * dist_matrix){
    _npop=npop, _ngen=ngen;
    
    int p1, p2; 
    // Read a pair of primes from file "Primes" to seed RNG
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    int seed[4]; // Array to store seed for RNG
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get current MPI process rank
    
    // Read seed values from "seed.in"
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    Seed.close();

    // Make the seed unique for each MPI process by adding rank-dependent offset
    for (int i = 0; i < 4; i++) {
        seed[i] += rank * 7; // 7 is an arbitrary constant chosen to differentiate RNG streams
    }

    // Initialize RNG with the unique seed and primes adjusted by rank to avoid correlations between ranks
    _rnd.SetRandom(seed, p1 + rank, p2 + rank);

    // Store the distance matrix by copying the provided pointer data
    _dist_matrix = *dist_matrix;

    _pop.set_size(_npop); // Resize population vector

    // Initialize each route in the population with the distance matrix and RNG
    for (int i = 0; i < _npop; ++i) {
        _pop[i].initialize(&_dist_matrix, _rnd);
    }
}

// Getter: return the Route at index i in the population
Route Population :: get_percorso(int i){
    return _pop(i);
}

// Select an index in the population based on a power-law distribution controlled by _selexp.
// This biases selection towards better individuals when _selexp > 1.
int Population::select_index() {
    double r = _rnd.Rannyu(); 
    int j = int(_npop * std::pow(r, _selexp)); 
    return j;
}

// Sort population by route length in ascending order, so best routes are at front
void Population::sort_by_length() {
    std::vector<Route> temp_vec(_npop);
    for (int i = 0; i < _npop; ++i){
        temp_vec[i] = _pop[i];
    }

    // Use std::sort with a lambda comparing route lengths by calculating with _dist_matrix
    std::sort(temp_vec.begin(), temp_vec.end(), [&](const Route& a, const Route& b) {
        return a.calculate_length(&_dist_matrix) < b.calculate_length(&_dist_matrix);
    });
    
    // Copy sorted routes back into _pop vector
    for (int i = 0; i < _npop; ++i){
        _pop[i] = temp_vec[i];
    }
}

// Select a route from population using select_index
Route Population::select() { 
    int j = select_index();
    return _pop(j);
}

// Return current population size
int Population::size() const {
    return _pop.size();
}

// Periodic boundary condition on city index; wraps index if >= ndim (110)
int Population :: pbc(int city) const{
    int ndim = 110;
    if(city >= ndim){
        return city - ndim + 1;
    }
    return city;
}

// Perform crossover between two parent routes to produce offspring routes
arma::field<Route> Population::crossover(Route a, Route b) {
    if (_rnd.Rannyu() < _pcross) { // Only crossover with probability _pcross
        int n_cities = a.getdim();
        arma::Col<int> parent1 = a.getroute();
        arma::Col<int> parent2 = b.getroute();

        // Choose two cut points for the crossover segment ensuring valid indices
        int cut1 = int(_rnd.Rannyu(1, n_cities - 2));
        int cut2 = int(_rnd.Rannyu(cut1 + 1, n_cities));

        // Initialize offspring vectors with zeros (assuming cities indexed > 0)
        arma::Col<int> offspring1(n_cities, arma::fill::zeros);
        arma::Col<int> offspring2(n_cities, arma::fill::zeros);

        // Copy segment between cut points from each parent to corresponding offspring
        for (int i = cut1; i <= cut2; ++i) {
            offspring1[i] = parent1[i];
            offspring2[i] = parent2[i];
        }

        // Fill the remaining positions in offspring1 with cities from parent2 in order, skipping duplicates
        int idx1 = (cut2 + 1) % n_cities;
        for (int i = 0; i < n_cities; ++i) {
            int city = parent2[(cut2 + 1 + i) % n_cities];
            if (arma::any(offspring1 == city)) continue; // Skip duplicates
            offspring1[idx1] = city;
            idx1 = (idx1 + 1) % n_cities;
        }

        // Similarly fill offspring2 with remaining cities from parent1
        int idx2 = idx1;
        for (int i = 0; i < n_cities; ++i) {
            int city = parent1[(cut2 + 1 + i) % n_cities];
            if (arma::any(offspring2 == city)) continue;
            offspring2[idx2] = city;
            idx2 = (idx2 + 1) % n_cities;
        }

        // Create Route objects for offspring and assign routes
        arma::field<Route> children(2);
        children[0].initialize(&_dist_matrix, _rnd);
        children[1].initialize(&_dist_matrix, _rnd);
        children[0].setroute(offspring1);
        children[1].setroute(offspring2);

        return children;
    } else {
        // No crossover: offspring are copies of parents
        return {a, b};
    }
}

// Sort vector a to match the order defined by reference vector ref
arma::Col<int> Population :: sort_by_reference(arma::Col<int> a, arma::Col<int> ref){
    arma::Col<int> ref_pos(ref.size()); 
    // Build a lookup table: position of each city in reference vector
    for(int i = 0; i < ref.size(); i++){
        ref_pos[ref[i]-1] = i;
    }

    // Convert a to std::vector for std::sort with custom comparator
    vector<int> sorted(a.begin(),a.end());
    sort(sorted.begin(), sorted.end(), [&ref_pos](int j, int k){
        return ref_pos[j] < ref_pos[k];
    });

    return arma::Col<int>(sorted); // Convert back to arma vector
}

// Main evolution function: run genetic algorithm for _ngen generations with migration every Nmigr generations
void Population::evolve(const arma::mat dist_matrix, int rank, int size, MPI_Comm comm, int Nmigr) {
    ofstream coutf;
    if(rank == 0) {
        // Only rank 0 writes results file with best and average route length per generation
        coutf.open("results.dat");
        coutf << "#" << "\t\t" << "L2" <<  "\t\t" << "<L2>" << endl;
    }

    sort_by_length(); // Sort population at start

    for (int i = 0; i < _ngen; i++) {
        arma::field<Route> figli(_npop);

        // Create offspring by selecting parents and performing crossover
        for (int j = 0; j < _npop / 2; j++) {
            Route a = select();
            Route b = select();
            arma::field<Route> temp = crossover(a, b);

            // Initialize offspring routes with distance matrix and RNG
            figli[j].initialize(&_dist_matrix, _rnd);
            figli[j + _npop / 2].initialize(&_dist_matrix, _rnd);

            // Set routes for offspring
            figli[j].setroute(temp[0].getroute());
            figli[j + _npop / 2].setroute(temp[1].getroute());
        }

        // Mutate all offspring routes
        for (int j = 0; j < _npop; j++) {
            figli[j].mutate();
        }

        // Replace old population with offspring
        _pop = figli;
        sort_by_length();

        // Every Nmigr generations, perform migration between MPI ranks
        if (i != 0 && i % Nmigr == 0) {
            migrate(rank, size, comm);
        }

        // Compute average length of best half of population
        double sum = 0;
        for (int j = 0; j < _npop / 2; j++) {
            sum += _pop[j].calculate_length(&dist_matrix);
        }

        if(rank == 0) {
            // Write current generation, best length, average length to file
            coutf << i << "\t\t" << _pop[0].calculate_length(&dist_matrix) << "\t\t" << sum / (_npop / 2) << endl;
        }
    }

    if(rank == 0) coutf.close();

    // Output the best route to a file (only rank 0)
    if(rank == 0) {
        ofstream out("best_route.dat");
        out << "#Best route \t \t x \t \t y" << endl;
        for (int i = 0; i < _pop[0].getdim(); i++) {
            out << _pop[0].getroute()[i] << endl;
        }
        out.close();
    }
}

// Replace the worst individual in the population with the migrant route, then re-sort population by length
void Population::replace_best(const Route& migrant) {
    _pop[_npop - 1] = migrant;
    sort_by_length();
}

// MPI migration function: exchange best routes between MPI processes
void Population :: migrate(int rank, int size, MPI_Comm comm) {
    Route best = _pop[0];  // Best route from this population
    arma::Col<int> send_data = best.getroute(); // Get route vector to send
    arma::Col<int> recv_data(send_data.n_elem); // Buffer for received route

    int ndim = send_data.n_elem;

    // Select random target process different from current rank
    int target;
    do {
        target = rand() % size;
    } while (target == rank);

    // Exchange best routes with the chosen target process via MPI_Sendrecv
    MPI_Sendrecv(
        send_data.memptr(), ndim, MPI_INT, target, 0,
        recv_data.memptr(), ndim, MPI_INT, target, 0,
        comm, MPI_STATUS_IGNORE
    );

    // Construct new route from received data
    Route migrant = best;
    migrant.setroute(recv_data);

    // Check validity of received route and replace worst if valid
    if (migrant.check()) {
        replace_best(migrant);
    } else {
        std::cerr << "[Rank " << rank << "] Received invalid route. Migration skipped.\n";
    }
}
