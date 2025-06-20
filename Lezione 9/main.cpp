#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include "popolazione.h"
#include "random.h"
#include "percorso.h"

using namespace std;

void generate_cities(int choice, int n_cities, arma::mat & dist_matrix);

int main(){

    int choice;

    cout << "Insert 0 for cities on a circumference, 1 for cities inside a square: ";
    cin >> choice;
    if (choice != 0 && choice != 1) {
        cout << "Invalid choice. Please enter 0 or 1." << endl;
        return 1; // Exit with error if invalid choice
    }
    
    int n_cities = 34;
    arma::mat dist_matrix(n_cities, n_cities);
    generate_cities(choice, n_cities, dist_matrix);

    Population pop;
    int npop = 200;   // Number of individuals in population
    int ngen = 1000;  // Number of generations for GA evolution
    pop.initialize_pop(npop, ngen, &dist_matrix); // Pass pointer to distance matrix
    pop.evolve(dist_matrix);

    return 0;
}

void generate_cities(int choice, int n_cities, arma::mat & dist_matrix) {
    ofstream outc("cities.dat");
    outc << "# \t \t x \t \t y" << endl;

    Random rnd_city;

    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2; // Read two prime numbers for RNG seeding
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd_city.SetRandom(seed, p1, p2); // Initialize RNG with seed and primes
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
    }

    vec x(n_cities);
    vec y(n_cities);

    if (choice == 0) {
        // Distribute cities evenly on circumference (unit circle)
        for (int i = 0; i < n_cities; i++) {
            double theta = rnd_city.Rannyu(0, 2 * M_PI);
            x[i] =  cos(theta);
            y[i] =  sin(theta);
            outc << i+1 << "\t\t" << x[i] << "\t\t" << y[i] << endl;
        }
    } else if (choice == 1) {
        // Distribute cities uniformly inside square [-1,1] x [-1,1]
        for (int i = 0; i < n_cities; i++) {
            x[i] = rnd_city.Rannyu(-1.0,1.0);
            y[i] = rnd_city.Rannyu(-1.0,1.0);
            outc << i+1 << "\t\t" << x(i) << "\t\t" << y(i) << endl;
        }
    }

    // Compute Euclidean distance matrix between all pairs of cities
    for(int i = 0; i < n_cities; i++) {
        for(int j = 0; j <  n_cities; j++) {
            dist_matrix(i,j) = sqrt(pow(x(i)-x(j),2) + pow(y(i)-y(j),2));
            // dist_matrix(j,i) = dist_matrix(i,j); // symmetric, but commented out
        }
    }

    outc.close();
}
