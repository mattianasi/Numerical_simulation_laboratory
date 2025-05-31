#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>
#include "random.h"

using namespace std;

double error(vector<double>, vector<double>, int);

int main() {
    Random rand;
    
    int N = 100;         // Number of blocks
    int M = 100000;      // Total number of simulations
    int L = M / N;       // Simulations per block
    double steps = 100;  // Number of steps in the random walk

    int seed[4];
    int p1, p2;

    // Initialize RNG from external files
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
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
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // 1. DISCRETE RANDOM WALK
    vector<double> sum_discrete(N, 0.0);
    vector<double> sum2_discrete(N, 0.0);
    vector<double> err_discrete(N, 0.0);
    double interval[6] = {1.0/6, 2.0/6, 3.0/6, 4.0/6, 5.0/6, 1.0}; // 6 directions

    for (int i = 0; i < N; i++) {
        vector<double> distance(steps, 0.0);
        for (int j = 0; j < L; j++) {
            vector<double> pos(6, 0.0);  // Discrete directions counter

            for (int k = 0; k < steps; k++) {
                double x = rand.Rannyu(); // Random number to determine direction
                for (int m = 0; m < 6; m++) {
                    if (x <= interval[m]) {
                        pos[m] += 1;
                        break;
                    }
                }
                // Distance squared from origin after k+1 steps
                distance[k] += pow(pos[0]-pos[1], 2) + pow(pos[2]-pos[3], 2) + pow(pos[4]-pos[5], 2);
            }
        }
        for (int k = 0; k < steps; k++) {
            distance[k] /= L;             // Average over L simulations
            distance[k] = sqrt(distance[k]); // Root mean square distance
            sum_discrete[k] += distance[k];
            sum2_discrete[k] += distance[k] * distance[k];
        }
    }

    ofstream out1("data.dat");
    for (int i = 0; i < steps; i++) {
        sum_discrete[i] /= N;
        sum2_discrete[i] /= N;
        err_discrete[i] = sqrt((sum2_discrete[i] - pow(sum_discrete[i], 2)) / N);
        out1 << sum_discrete[i] << "," << err_discrete[i] << endl;
    }
    out1.close();

    // 2. CONTINUOUS RANDOM WALK
    vector<double> sum_cont(N, 0.0);
    vector<double> sum2_cont(N, 0.0);
    vector<double> err_cont(N, 0.0);

    for (int i = 0; i < N; i++) {
        vector<double> distance(steps, 0.0);

        for (int j = 0; j < L; j++) {
            vector<double> r(3, 0.0); // x, y, z

            for (int k = 0; k < steps; k++) {
                double a = rand.Rannyu();
                double phi = rand.Rannyu(0, 2*M_PI);
                double theta = acos(1 - 2*a); // Uniform over sphere

                // Step in spherical coordinates
                r[0] += sin(theta) * cos(phi);
                r[1] += sin(theta) * sin(phi);
                r[2] += cos(theta);

                distance[k] += r[0]*r[0] + r[1]*r[1] + r[2]*r[2]; // r²
            }
        }

        for (int k = 0; k < steps; k++) {
            distance[k] /= L;
            distance[k] = sqrt(distance[k]);
            sum_cont[k] += distance[k];
            sum2_cont[k] += distance[k] * distance[k];
        }
    }

    ofstream out2("datacont.dat");
    for (int i = 0; i < steps; i++) {
        sum_cont[i] /= N;
        sum2_cont[i] /= N;
        err_cont[i] = sqrt((sum2_cont[i] - pow(sum_cont[i], 2)) / N);
        out2 << sum_cont[i] << "," << err_cont[i] << endl;
    }
    out2.close();

    return 0;
}

// Standard statistical error calculation
double error(vector<double> AV, vector<double> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
    }
}
