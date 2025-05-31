#include "funzioni.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

// Function to load primes from a file
void LoadPrimes(int &p1, int &p2) {
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
        Primes.close();
    } else {
        cerr << "PROBLEM: Unable to open Primes file" << endl;
        exit(1);
    }
}

// Function to load the seed for random number generation
void LoadSeed(int seed[4], Random &rand, int p1, int p2) {
    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rand.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        exit(1);
    }
}

// Function to calculate the error for the averages
double error(const vector<double>& averages, const vector<double>& averages2, int i) {
    if (i == 0) return 0.0;
    return sqrt((averages2[i] - pow(averages[i], 2)) / i);
}

// Function to calculate the distance between points
double CalculateDistance(const vector<int>& pos) {
    return sqrt(pow((pos[0] - pos[1]), 2) + pow((pos[2] - pos[3]), 2) + pow((pos[4] - pos[5]), 2));
}