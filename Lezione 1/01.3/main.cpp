#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"

using namespace std;

// Function declarations
double error(vector<double> AV, vector<double> AV2, int n);
bool Intersect(double l, double d, double ang, Random& rnd);

int main() {
    Random rand;
    int seed[4];
    int p1, p2;

    // Reading prime numbers for random number generator initialization
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    // Reading the initial seed
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

    // ~~~~ Estimation of pi using Buffon's Needle ~~~~ //

    int M = 1000000; // Total number of throws
    int N = 100;     // Number of blocks
    int L = M / N;   // Number of throws per block
    int l = 1;       // Needle length
    double d = 3;    // Distance between lines

    vector<double> sum_pi(N), sum_pi2(N), err_pi(N);
    double theta;

    // Main simulation loop
    for (int i = 0; i < N; i++) {
        int Nhit = 0; // Number of intersections in this block
        for (int j = 0; j < L; j++) {
            // Generate a random direction using polar coordinates
            double x = rand.Rannyu(-1, 1);
            double y = rand.Rannyu(-1, 1);
            if (y >= 0) {
                theta = acos(x / sqrt(x * x + y * y));
            } else {
                theta = 2 * M_PI - acos(x / sqrt(x * x + y * y));
            }

            if (Intersect(l, d, theta / 2, rand)) {
                Nhit++;
            }
        }
        sum_pi[i] = (2.0 * l * L) / (d * Nhit);   // Block estimate of pi
        sum_pi2[i] = pow(sum_pi[i], 2);          // Squared value for variance estimation
    }

    // Output results
    ofstream out("pi.dat");
    vector<double> sum_prog(N, 0.0), sum_prog2(N, 0.0);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += sum_pi[j];
            sum_prog2[i] += sum_pi2[j];
        }
        sum_prog[i] /= (i + 1);
        sum_prog2[i] /= (i + 1);
        err_pi[i] = error(sum_prog, sum_prog2, i);

        out << sum_prog[i] << "\t" << err_pi[i] << endl;
    }
    out.close();

    cout << "end program" << endl;
    return 0;
}

// Function to compute statistical uncertainty
double error(vector<double> AV, vector<double> AV2, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
    }
}

// Function to determine whether the needle intersects a line
bool Intersect(double l, double d, double ang, Random& rnd) {
    double y = (l / 2) * sin(ang);      // Vertical projection of the half-needle
    double x = rnd.Rannyu() * d;        // Distance from center to nearest line
    return y >= min(x, d - x);          // True if the needle crosses a line
}
