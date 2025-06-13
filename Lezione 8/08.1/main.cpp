#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include <algorithm>
#include <filesystem>
#include "random.h"

using namespace std;

double error( double AV, double AV2 , int n );
double psi(double x, double mu, double sigma);
bool metro(double x , double y, double mu, double sigma, Random&  rnd);
double Hpsi(double x , double mu, double sigma);
double set_delta(Random& rnd, double& delta, double mu, double sigma, double x0, int steps);


int main (int argc, char *argv[]){

    // random number generator
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
 
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;


    // Simulation parameters
    double x = 0.0; // initial configuration (position)
    double mu = 0.772229;
    double sigma = 0.618926;
    	 	
    double delta = 0.5; // start amplitude of Metropolis step

    int n_steps = 10000; // number of steps per block
    int n_blocks = 100; // number of blocks
    int tune_steps = 1000; // number of steps to tune delta

    // Tune delta to reach target acceptance
    double tuned_acceptance = set_delta(rnd, delta, mu, sigma, x, tune_steps);

    ofstream coute("energy.dat");
    ofstream couta("acceptance.dat");
    ofstream amp("amplitude.dat");
    coute << "#     BLOCK:  ACTUAL_E:     E_AVE:      ERROR:" << endl;
    couta << "# BLOCK\tACCEPTANCE\n";
    amp << "# x\t|psi(x)|^2\n";
    
    double av; // holds the progressive average of the first (i+1) blocks
    double av2; // holds the progressive average of the squares of the first (i+1) blocks
    double sum_prog = 0.; // accumulates sum of averages of first (i+1) blocks
    double sum_prog2 = 0.; // accumulates sum of squared averages 

    double n_attemps;
    double n_accepted;

    for(int i=0; i <n_blocks; i++){ // loop over blocks
        double stima = 0;
        n_attemps = 0;
        n_accepted = 0;
        for(int j=0; j < n_steps; j++){ // loop over steps in a block
            // Metropolis step 
            n_attemps++;
            double y = x +  rnd.Rannyu(-1.0,1.0) * delta; // propose new position
            if(metro(x,y,mu, sigma,rnd)){ // accept or reject move
                x = y;
                n_accepted++;
            }
            stima +=  Hpsi(x, mu, sigma); // accumulate local energy estimate
            amp << x << "\t" << pow(psi(x, mu, sigma), 2) << endl; // record |psi(x)|^2
        }
        couta << (i+1) << "\t" << double(n_accepted)/double(n_attemps) << endl; // acceptance per block
        sum_prog= sum_prog + stima/n_steps;
        sum_prog2 =sum_prog2 + stima/n_steps * stima/n_steps;
        av = (sum_prog) / (i+1);
        av2 = ( sum_prog2 ) / (i+1);
        coute << (i+1) << "\t" << stima/n_steps  << "\t" << (av ) << "\t" << error(av, av2, i) << endl;
        }
    coute.close();
    couta.close();      

  return 0;
}


double error( double AV, double AV2 , int n ){
    // function returning the standard deviation of n elements
    if(n == 0){
        return 0;
    }
    // standard error of the mean
    return sqrt( (AV2 - AV*AV)/n );
    // note: usually divided by (n-1), but here divided by n since index starts at 0
}

double psi(double x, double mu, double sigma){
    // Trial wavefunction: sum of two Gaussians centered at ±mu
    return exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
}

double Hpsi(double x , double mu, double sigma){
    // Local energy = (H psi)/psi = kinetic + potential parts
    return -0.5 *(  
        -1/(sigma*sigma) *  exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) 
        -1/(sigma*sigma) *  exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) 
        + 1/(pow(sigma,4)) * (x-mu)*(x-mu) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) 
        + 1/(pow(sigma,4)) * (x+mu)*(x+mu) * exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) 
        ) / psi(x, mu, sigma) 
        + pow(x,4) - 5./2. * pow(x,2);
}

bool metro(double x , double y, double mu, double sigma, Random&  rnd){ // Metropolis algorithm
    bool decision = false;
    // acceptance probability proportional to squared wavefunction ratio
    double acceptance = min( 1. , psi(y , mu , sigma)*psi(y , mu , sigma) / (psi(x, mu , sigma)*psi(x, mu,sigma)) ); 
    if(rnd.Rannyu() < acceptance ) decision = true; // accept move with Metropolis criterion
    return decision;
}

double set_delta(Random& rnd, double& delta, double mu, double sigma, double x0, int steps) {
    double target_acceptance=0.5; // target acceptance ratio
    double acceptance = 0;
    double x = x0; // starting position
    double tol = 0.1; // tolerance for acceptance to be close to target
    int max_iter = 1000;
    int iter = 0;
    while (iter < max_iter) {
        int accepted = 0;
        for (int i = 0; i < steps; ++i) {
            double y = x + rnd.Rannyu(-1., 1.) * delta; // propose step
            if (metro(x, y, mu, sigma, rnd)) {
                x = y;
                accepted++;
            }
        }
        acceptance = static_cast<double>(accepted) / steps;
        if (fabs(acceptance - target_acceptance) < tol) {
            break; // exit if acceptance is close enough to target
        }
        if (acceptance < target_acceptance)
            delta *= 0.9; // reduce delta if acceptance too low
        else
            delta *= 1.1; // increase delta if acceptance too high

        iter++;
    }
    return acceptance;
}
