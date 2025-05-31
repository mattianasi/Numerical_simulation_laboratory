#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"

using namespace std;

int main() {
    Random rand;
    int N = 10000;

    // Vectors to store the averages of samples
    vector<double> S1(N);     // Single sample
    vector<double> S2(N);     // Average of 2 samples
    vector<double> S10(N);    // Average of 10 samples
    vector<double> S100(N);   // Average of 100 samples

    ofstream out1, out2, out3;  // Output files

    // Initialize random number generator with seed
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
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
                rand.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
    }

    // Uniform distribution
    out1.open("lin.dat");
    for (int i = 0; i < N; i++) {
        S10[i] = 0;
        S100[i] = 0;

        S1[i] = rand.Rannyu();                             // 1 sample
        S2[i] = (rand.Rannyu() + rand.Rannyu()) / 2;       // 2 samples average
        for (int j = 0; j < 10; j++) {
            S10[i] += rand.Rannyu();                       // Sum 10 samples
        }
        for (int j = 0; j < 100; j++) {
            S100[i] += rand.Rannyu();                      // Sum 100 samples
        }

        // Write normalized averages
        out1 << S1[i] << "," << S2[i] << "," << S10[i] / 10 << "," << S100[i] / 100 << endl;
    }
    out1.close();

    // Exponential distribution
    out2.open("exp.dat");
    for (int i = 0; i < N; i++) {
        S10[i] = 0;
        S100[i] = 0;

        S1[i] = rand.Exp(1);                              // 1 sample
        S2[i] = (rand.Exp(1) + rand.Exp(1)) / 2;          // 2 samples average
        for (int j = 0; j < 10; j++) {
            S10[i] += rand.Exp(1);                        // Sum 10 samples
        }
        for (int j = 0; j < 100; j++) {
            S100[i] += rand.Exp(1);                       // Sum 100 samples
        }

        out2 << S1[i] << "," << S2[i] << "," << S10[i] / 10 << "," << S100[i] / 100 << endl;
    }
    out2.close();

    // Lorentzian (Cauchy) distribution
    out3.open("lorentz.dat");
    for (int i = 0; i < N; i++) {
        S10[i] = 0;
        S100[i] = 0;

        S1[i] = rand.Lorentz(0, 1);                             // 1 sample
        S2[i] = (rand.Lorentz(0, 1) + rand.Lorentz(0, 1)) / 2;  // 2 samples average
        for (int j = 0; j < 10; j++) {
            S10[i] += rand.Lorentz(0, 1);                       // Sum 10 samples
        }
        for (int j = 0; j < 100; j++) {
            S100[i] += rand.Lorentz(0, 1);                      // Sum 100 samples
        }

        out3 << S1[i] << "," << S2[i] << "," << S10[i] / 10 << "," << S100[i] / 100 << endl;
    }
    out3.close();

    return 0;
}
