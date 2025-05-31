#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "random.h"

using namespace std;

// Function to calculate statistical error of the mean
double error(vector<double>, vector<double>, int);

// Functions to integrate
double function(double);
double function1(double);

int main(){

    Random rand;  // Random number generator object

    int seed[4];
    int p1, p2;
    // Read prime numbers needed to initialize RNG
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    // Read seed values from file and initialize RNG
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rand.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    int M = 100000; // Total number of random samples
    int N = 100;    // Number of blocks for blocking method
    int L = M/N;    // Number of samples per block

    vector<double> sum(N);   // Block averages
    vector<double> sum2(N);  // Square of block averages for variance calculation
    vector<double> err(N);   // Errors from blocking method
    ofstream out1,out2;

    // First integration: using uniform random numbers (rand.Rannyu())
    for(int i=0; i<N; i++){
        double f = 0;
        for(int j=0; j<L; j++){
            double x = rand.Rannyu();   // Uniform random number in [0,1)
            f += function(x);           // Evaluate integrand at x
        }
        sum[i] = f/L;                   // Average value of function in block i
        sum2[i] = sum[i]*sum[i];       // Store square for error estimation
    }

    // Open output file for first integral
    out1.open("data.dat");

    vector<double> averages(N, 0.0);   // Progressive averages over blocks
    vector<double> averages2(N, 0.0);  // Progressive averages of squared values

    // Calculate progressive averages and errors using blocking method
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i + 1; j++) {
            averages[i] += sum[j];
            averages2[i] += sum2[j];
        }
        averages[i] /= (i + 1);
        averages2[i] /= (i + 1);
        err[i] = error(averages, averages2, i);
        out1 << averages[i] << "," << err[i] << endl;  // Write average and error
    }
    out1.close();

    // Second integration: using a different distribution (rand.Distr())
    vector<double> sum_(N);
    vector<double> sum2_(N);
    vector<double> err_(N);

    for(int i=0; i<N; i++){
        double f = 0;
        for(int j=0; j<L; j++){
            double x = rand.Distr();      // Random number from custom distribution
            f += function1(x);            // Evaluate integrand adapted for this distribution
        }
        sum_[i] = f/L;
        sum2_[i] = sum_[i]*sum_[i];
    }

    // Output file for second integral
    out2.open("data1.dat");

    vector<double> averages_(N, 0.0);
    vector<double> averages2_(N, 0.0);

    // Calculate progressive averages and errors for second integral
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i + 1; j++) {
            averages_[i] += sum_[j];
            averages2_[i] += sum2_[j];
        }
        averages_[i] /= (i + 1);
        averages2_[i] /= (i + 1);
        err_[i] = error(averages_, averages2_, i);
        out2 << averages_[i] << "," << err_[i] << endl;
    }
    out2.close();

    return 0;
}

// Function to compute statistical error for progressive averages
double error(vector<double> AV, vector<double> AV2, int n){
    if(n == 0) {
        return 0;  // No error for first point
    } else {
        // Standard error of the mean formula
        return sqrt((AV2[n]-pow(AV[n],2))/n);
    }
}

// Function to be integrated in first method
double function(double x){
    return M_PI/2*cos(M_PI*x/2);
}

// Function to be integrated in second method, designed to work with importance sampling, calculated from Taylor expansion 
double function1(double x){
    return M_PI/4*(cos(M_PI*x/2))/(1-x);
}
