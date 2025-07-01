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

void load_cities(mat& dist_matrix);

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int ncap = 110;
    const int npop = 100;
    const int ngen = 100000;
    const int Nmigr = 500; // Migration frequency

    mat dist_matrix(ncap, ncap);
    load_cities(dist_matrix);

    Population pop;
    pop.initialize_pop(npop, ngen, &dist_matrix);

    // Pass rank, size, communicator and migration frequency to evolve
    pop.evolve(dist_matrix, rank, size, MPI_COMM_WORLD, Nmigr);

    // 1) Get local best length and best route
    double local_best_length = pop.get_best_length();
    arma::Col<int> local_best_route = pop.get_best_route();

    int n_cities = local_best_route.n_elem;

    // 2) Gather all local best lengths to rank 0
    double* all_lengths = nullptr;
    if (rank == 0) {
        all_lengths = new double[size];
    }

    MPI_Gather(&local_best_length, 1, MPI_DOUBLE,
            all_lengths, 1, MPI_DOUBLE,
            0, MPI_COMM_WORLD);

    // 3) Rank 0 finds the global best and which rank owns it
    int best_owner_rank = 0;
    if (rank == 0) {
        double min_length = all_lengths[0];
        for (int i = 1; i < size; i++) {
            if (all_lengths[i] < min_length) {
                min_length = all_lengths[i];
                best_owner_rank = i;
            }
        }
        delete[] all_lengths;
    }

    // 4) Broadcast best_owner_rank to all
    MPI_Bcast(&best_owner_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 5) Send best route from best_owner_rank to rank 0 (if best_owner_rank != 0)
    if (rank == best_owner_rank && best_owner_rank != 0) {
        MPI_Send(local_best_route.memptr(), n_cities, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    arma::Col<int> best_route_global(n_cities);
    if (rank == 0) {
        if (best_owner_rank == 0) {
            best_route_global = local_best_route;
        } else {
            MPI_Recv(best_route_global.memptr(), n_cities, MPI_INT, best_owner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Write the best route to file
        ofstream out("best_route.dat");
        out << "#Best route \t \t x \t \t y" << endl;
        for (int i = 0; i < n_cities; i++) {
            out << best_route_global[i] << endl;
        }
        out.close();
    }


    MPI_Finalize();
    return 0;
}

void load_cities(mat& dist_matrix) {
    int ncap = 110;
    vec x(ncap), y(ncap);
    ifstream filein("cap_prov_ita.dat");
    if (!filein) {
        cerr << "Error: could not open cap_prov_ita.dat" << endl;
        exit(1);
    }
    for (int i = 0; i < ncap; i++) {
        filein >> x(i) >> y(i);
    }
    filein.close();

    for (int i = 0; i < ncap; i++) {
        for (int j = 0; j < ncap; j++) {
            dist_matrix(i,j) = pow(x(i)-x(j),2) + pow(y(i)-y(j),2);
        }
    }
}
