#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "random.h" 

using namespace std;

// Function that calculates the statistical error
double error(vector<double> AV, vector<double> AV2, int n);

int main() {
    // Simulation parameters
    int M = 1000000;  // Total number of random throws
    int N = 100;      // Number of blocks
    int L = M / N;    // Number of throws per block

    Random rand;      // Object for random number generation
    vector<double> av, av2; // Vectors for block averages and squares of block averages
    vector<double> sum_prog(N), sum2_prog(N), err_prog(N); // Progressive sums and errors
    ofstream out1, out2, out3; // Output files for results

    // Read the seed from the file containing prime numbers
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2; // Read two prime numbers
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl; // Error message if the file can't be opened
    }
    Primes.close();

    // Read the random number generator seed from seed.in
    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rand.SetRandom(seed, p1, p2); // Initialize the RNG
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl; // Error message if the file can't be opened
    }

    // Compute the average of L random numbers and its variance (Part 1)
    for (int i = 0; i < N; i++) {
        double sum1 = 0;
        for (int j = 0; j < L; j++) {
            sum1 += rand.Rannyu(); // Sum the random numbers
        }
        av.push_back(sum1 / L);   // Block average
        av2.push_back(pow(av[i], 2)); // Square of the block average
    }

    // Write block averages and their errors to output files
    out1.open("media.dat");
    out2.open("errore.dat");

    for (int i = 0; i < N; i++) {
        // Progressive sum for average and its square
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += av[j];
            sum2_prog[i] += av2[j];
        }
        sum_prog[i] /= (i + 1); // Progressive average
        out1 << sum_prog[i] << endl;
        sum2_prog[i] /= (i + 1); // Progressive average of squares
        err_prog[i] = error(sum_prog, sum2_prog, i); // Calculate error
        out2 << err_prog[i] << endl;
    }

    out1.close();
    out2.close();

    // Compute variance and its error (Part 2)
    vector<double> av_, av2_; // Vectors for Part 2
    vector<double> sum_prog_(N), sum2_prog_(N), err_prog_(N);

    for (int i = 0; i < N; i++) {
        double sum2 = 0;
        for (int j = 0; j < L; j++) {
            sum2 += pow(rand.Rannyu() - 0.5, 2); // Variance of uniform distribution in [0,1]
        }
        av_.push_back(sum2 / L);   // Block average of variance
        av2_.push_back(pow(av_[i], 2)); // Square of the block average
    }

    out1.open("var.dat");
    out2.open("errorevar.dat");

    for (int i = 0; i < N; i++) {
        // Progressive sum for average and its square
        for (int j = 0; j <= i; j++) {
            sum_prog_[i] += av_[j];
            sum2_prog_[i] += av2_[j];
        }
        sum_prog_[i] /= (i + 1); // Progressive average
        out1 << sum_prog_[i] << endl;
        sum2_prog_[i] /= (i + 1); // Progressive average of squares
        err_prog_[i] = error(sum_prog_, sum2_prog_, i); // Calculate error
        out2 << err_prog_[i] << endl;
    }

    out1.close();
    out2.close();

    // Compute Chi-squared values for uniform distribution test
    int m = 100; // Number of intervals
    int n = 10000; // Number of throws
    double expected = double(n) / m; // Expected count in each interval
    vector<double> chi2(m); // Vector to store Chi-squared values
    out3.open("chi2.dat");

    for (int i = 0; i < m; i++) {
        vector<double> ni(m, 0); // Vector for observed frequencies
        double chi2_ = 0;

        for (int j = 0; j < n; j++) {
            double x = rand.Rannyu(); // Generate a random number
            int interval = floor(x * m); // Determine which interval it falls into
            ni[interval]++; // Increment observed frequency
        }

        // Compute Chi-squared
        for (int j = 0; j < m; j++) {
            chi2_ += pow(ni[j] - expected, 2) / expected;
        }

        chi2[i] = chi2_; // Save Chi-squared value
        out3 << chi2_ << endl;
    }

    out3.close();

    // Save the seed to reproduce the simulation
    rand.SaveSeed();

    return 0;
}

// Function to compute statistical error (standard deviation)
double error(vector<double> AV, vector<double> AV2, int n) {
    if (n == 0) {
        return 0; // Error is 0 at the first step
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2)) / n); // Compute the error
    }
}
