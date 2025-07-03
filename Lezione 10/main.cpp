#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include "popolazione.h"
#include "random.h"
#include "percorso.h"
#include <armadillo>

using namespace std;
using namespace arma;

// Function declaration to load city coordinates and fill distance matrix
void load_cities(mat& dist_matrix);

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);  // Initialize MPI environment

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get MPI process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Get total number of MPI processes

    // Problem and GA parameters
    const int ncap = 110;     // Number of cities
    const int npop = 100;     // Population size
    const int ngen = 10000;   // Number of generations to evolve
    const int Nmigr = 50;     // Frequency of migration between MPI processes

    mat dist_matrix(ncap, ncap);  // Matrix storing distances between cities
    load_cities(dist_matrix);      // Load city coordinates and compute distance matrix

    Population pop;                      // Create population object
    pop.initialize_pop(npop, ngen, &dist_matrix);  // Initialize population with given parameters

    // Run evolution process: pass MPI rank, size, communicator and migration frequency
    pop.evolve(dist_matrix, rank, size, MPI_COMM_WORLD, Nmigr);

    // After evolution, get local best route length and the corresponding route for this MPI process
    double local_best_length = pop.get_best_length();
    arma::Col<int> local_best_route = pop.get_best_route();

    int n_cities = local_best_route.n_elem;  // Number of cities in the route

    // Prepare to gather all local best lengths at rank 0
    double* all_lengths = nullptr;
    if (rank == 0) {
        all_lengths = new double[size];  // Allocate array to hold lengths from all processes
    }

    // Gather best lengths from all processes to rank 0
    MPI_Gather(&local_best_length, 1, MPI_DOUBLE,
               all_lengths, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // Rank 0 determines which process has the globally best (shortest) route
    int best_owner_rank = 0;
    if (rank == 0) {
        double min_length = all_lengths[0];
        for (int i = 1; i < size; i++) {
            if (all_lengths[i] < min_length) {
                min_length = all_lengths[i];
                best_owner_rank = i;
            }
        }
        delete[] all_lengths;  // Clean up allocated memory
    }

    // Broadcast the rank owning the best route to all MPI processes
    MPI_Bcast(&best_owner_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // If current process owns the best route but is not rank 0, send the best route to rank 0
    if (rank == best_owner_rank && best_owner_rank != 0) {
        MPI_Send(local_best_route.memptr(), n_cities, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    arma::Col<int> best_route_global(n_cities);  // Container for global best route on rank 0

    if (rank == 0) {
        if (best_owner_rank == 0) {
            best_route_global = local_best_route;  // Best route is local to rank 0
        } else {
            // Receive best route from best_owner_rank process
            MPI_Recv(best_route_global.memptr(), n_cities, MPI_INT, best_owner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Write the global best route to a file for later inspection/plotting
        ofstream out("best_route.dat");
        out << "#Best route \t \t x \t \t y" << endl;
        for (int i = 0; i < n_cities; i++) {
            out << best_route_global[i] << endl;
        }
        out.close();
    }

    MPI_Finalize();  // Finalize MPI environment before exit
    return 0;
}

// Load city coordinates from file and compute squared Euclidean distance matrix
void load_cities(mat& dist_matrix) {
    int ncap = 110;
    vec x(ncap), y(ncap);

    ifstream filein("cap_prov_ita.dat");  // Open file with city coordinates
    if (!filein) {
        cerr << "Error: could not open cap_prov_ita.dat" << endl;
        exit(1);  // Exit if file not found
    }

    // Read x and y coordinates for all cities
    for (int i = 0; i < ncap; i++) {
        filein >> x(i) >> y(i);
    }
    filein.close();

    // Compute squared Euclidean distances between every pair of cities
    for (int i = 0; i < ncap; i++) {
        for (int j = 0; j < ncap; j++) {
            dist_matrix(i,j) = pow(x(i) - x(j), 2) + pow(y(i) - y(j), 2);
        }
    }
}
