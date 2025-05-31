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
        return 1; // Exit the program with an error code
    }
    
    int n_cities = 34; // Number of cities
    arma::mat dist_matrix(n_cities, n_cities); // Distance matrix
    generate_cities(choice, n_cities, dist_matrix); // Generate cities and distance matrix

    Population pop; // Create a population object
    int npop = 200; // Number of individuals in the population
    int ngen = 1000; // Number of generations
    pop.initialize_pop(npop, ngen, &dist_matrix); // Initialize the population
    pop.evolve(dist_matrix); // Evolve the population


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
        Primes >> p1 >> p2; // Legge due numeri primi
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl; // Messaggio di errore se il file non si apre
    }
    Primes.close();

    // Leggi il seme per la generazione di numeri casuali dal file seed.in
    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd_city.SetRandom(seed, p1, p2); // Inizializza il generatore di numeri casuali
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl; // Messaggio di errore se il file non si apre
    }

    vec x(n_cities);
    vec y(n_cities);

    if (choice == 0) {
        // Generate cities on a circumference
        for (int i = 0; i < n_cities; i++) {
            double theta = rnd_city.Rannyu(0, 2 * M_PI);
            x[i] =  cos(theta);
            y[i] =  sin(theta);
            outc << i+1 << "\t\t" << x[i] << "\t\t" << y[i] << endl;
        }
    } else if (choice == 1) {
        // Generate cities inside a square
        for (int i = 0; i < n_cities; i++) {
            x[i] = rnd_city.Rannyu(-1.0,1.0); // Random x coordinate in [-1, 1]
            y[i] = rnd_city.Rannyu(-1.0,1.0); // Random y coordinate in [-1, 1]
            outc << i+1 << "\t\t" << x(i) << "\t\t" << y(i) << endl;
        }
    }

    // Calculate the distance matrix
    for(int i = 0; i < n_cities; i++) {
        for(int j = 0; j <  n_cities; j++) {
            dist_matrix(i,j) = sqrt(pow(x(i)-x(j),2) + pow(y(i)-y(j),2)); // L1 norm
            //dist_matrix(j,i) = dist_matrix(i,j); // Symmetric matrix
        }
    }

    outc.close();
}



