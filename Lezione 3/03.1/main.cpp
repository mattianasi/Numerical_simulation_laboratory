#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"    // Custom random number generator class

using namespace std;

// Function to calculate the statistical error (standard error of the mean)
double error(vector<double>, vector<double>, int);

int main(){
    Random rand;
    int seed[4];
    int p1, p2;

    // Read two prime numbers from the file "Primes" to initialize the RNG
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    // Read seed values from the file "seed.in" and initialize the RNG
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while (!input.eof()){
            input >> property;
            if (property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rand.SetRandom(seed, p1, p2);  // Set RNG seed and primes
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // Parameters for Monte Carlo simulation of option pricing
    int M = 100000;   // total number of asset prices
    int N = 100;      // number of blocks for statistics
    int L = M / N;    // number of simulations per block

    // Option and market parameters
    double S0 = 100;    // initial stock price
    double T = 1;       // maturity (1 year)
    double K = 100;     // strike price
    double r = 0.1;     // risk-free interest rate
    double sigma = 0.25;// volatility

    // Vectors to store averages, squared averages, and errors
    vector<double> av(N), av2(N), err(N);

    // --------- First part: European Call Option (direct sampling) ---------
    for(int i = 0; i < N; i++){
        double C = 0;
        for(int j = 0; j < L; j++){
            double z = rand.Gauss(0,1);  // sample standard normal variable
            // Simulate stock price at maturity using geometric Brownian motion formula
            double S = S0 * exp((r - sigma*sigma/2)*T + sigma*z*sqrt(T));
            // Payoff discounted to present value (max between S-K and 0)
            C += exp(-r*T) * fmax(0, S - K);
        }
        av[i] = C / L;         // Average payoff for block i
        av2[i] = pow(av[i], 2); // Square of the average (for error calculation)
    }

    // Open output file to write results of call option
    ofstream out("data1.dat");
    vector<double> sum_prog(N, 0.0), sum_prog2(N, 0.0);

    // Calculate progressive averages and errors for each block
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += av[j];    // sum of averages up to block i
            sum_prog2[i] += av2[j];  // sum of squared averages up to block i
        }
        sum_prog[i] /= (i + 1);       // progressive average
        sum_prog2[i] /= (i + 1);      // progressive squared average
        err[i] = error(sum_prog, sum_prog2, i);  // calculate statistical error
        out << sum_prog[i] << "\t" << err[i] << endl;  // write results to file
    }
    out.close();

    // --------- Second part: European Put Option (direct sampling) ---------
    for(int i = 0; i < N; i++){
        double P = 0;
        for(int j = 0; j < L; j++){
            double z = rand.Gauss(0,1);
            double S = S0 * exp((r - sigma*sigma/2)*T + sigma*z*sqrt(T));
            P += exp(-r*T) * fmax(0, K - S);
        }
        av[i] = P / L;
        av2[i] = pow(av[i], 2);
    }

    out.open("data2.dat");

    // Calculate progressive averages and errors again
    for (int i = 0; i < N; i++) {
        sum_prog[i] = 0;
        sum_prog2[i] = 0;
        err[i] = 0;
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += av[j];
            sum_prog2[i] += av2[j];
        }
        sum_prog[i] /= (i + 1);
        sum_prog2[i] /= (i + 1);
        err[i] = error(sum_prog, sum_prog2, i);
        out << sum_prog[i] << "\t" << err[i] << endl;
    }
    out.close();

    // --------- Third part: European Call Option (discrete sampling) ---------
    // This simulates the stock price evolution in 100 discrete time steps
    for(int i = 0; i < N; i++){
        double C = 0;
        for(int j = 0; j < L; j++){
            double S1 = S0;
            double S2;
            for(int k = 0; k < 100; k++){
                double z = rand.Gauss(0,1);
                // Simulate price at next time step
                S2 = S1 * exp((r - sigma*sigma/2) * (1.0/100) + sigma * z * sqrt(1.0/100));
                S1 = S2;
            }
            // Payoff at maturity discounted to present
            C += exp(-r*T) * fmax(0, S2 - K);
        }
        av[i] = C / L;
        av2[i] = pow(av[i], 2);
    }

    out.open("data3.dat");

    // Progressive averages and errors
    for (int i = 0; i < N; i++) {
        sum_prog[i] = 0;
        sum_prog2[i] = 0;
        err[i] = 0;
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += av[j];
            sum_prog2[i] += av2[j];
        }
        sum_prog[i] /= (i + 1);
        sum_prog2[i] /= (i + 1);
        err[i] = error(sum_prog, sum_prog2, i);
        out << sum_prog[i] << "\t" << err[i] << endl;
    }
    out.close();

    // --------- Fourth part: European Put Option (discrete sampling) ---------
    for(int i = 0; i < N; i++){
        double P = 0;
        for(int j = 0; j < L; j++){
            double S1 = S0;
            double S2;
            for(int k = 0; k < 100; k++){
                double z = rand.Gauss(0,1);
                S2 = S1 * exp((r - sigma*sigma/2) * (1.0/100) + sigma * z * sqrt(1.0/100));
                S1 = S2;
            }
            P += exp(-r*T) * fmax(0, K - S2);
        }
        av[i] = P / L;
        av2[i] = pow(av[i], 2);
    }

    out.open("data4.dat");

    // Progressive averages and errors
    for (int i = 0; i < N; i++) {
        sum_prog[i] = 0;
        sum_prog2[i] = 0;
        err[i] = 0;
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += av[j];
            sum_prog2[i] += av2[j];
        }
        sum_prog[i] /= (i + 1);
        sum_prog2[i] /= (i + 1);
        err[i] = error(sum_prog, sum_prog2, i);
        out << sum_prog[i] << "\t" << err[i] << endl;
    }
    out.close();

    cout << "end program" << endl;
    return 0;
}

// Function to calculate the statistical error (standard error of the mean)
// AV = vector of progressive averages
// AV2 = vector of progressive squared averages
// n = index of the current block
double error(vector<double> AV, vector<double> AV2, int n){
    if(n == 0) {
        return 0;  // no error for first data point
    }
    else{
        // Standard error formula: sqrt((<x^2> - <x>^2) / n)
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
    }
};
