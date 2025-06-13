#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include <algorithm>
#include <vector>
#include "random.h"

using namespace std;

double error( double AV, double AV2 , int n );
double psi(double x, double mu, double sigma);
bool metro_H(double x , double y, double mu, double sigma, Random&  rnd);
double Hpsi(double x , double mu, double sigma);
double set_delta(Random& rnd, double& delta, double mu, double sigma, double x0, int steps);
vector<double> compute_H ( Random& rnd , double mu, double sigma );
bool metro_par(Random& rnd, double H_old , double H_new , double beta );


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

    double mu = 1.0, sigma = 1.0;
    vector <double> H = compute_H( rnd, mu, sigma); // initial parameters

    ofstream coutp("parameters.dat");
    coutp << "T  \t \t mu \t \t sigma \t \t H \t \t err" << endl;

    for ( double T = 2. ; T>= 0.01  ; T*=0.99){ // loop over temperature decreasing
        for ( int n = 0 ; n < 100 ; n++){
            double delta_mu = 0.5 * T; // Metropolis step amplitude (modified based on temperature)
            double delta_sigma = 0.5 * T; // Metropolis step amplitude (modified based on temperature)
            double beta = 1./T; // beta
            // do-while to control sigma, skip too low values
            double mu_new = fabs(mu + rnd.Rannyu(-1, 1)*delta_mu); // new proposal for mu
            double sigma_new;
            do{
                sigma_new = fabs(sigma + rnd.Rannyu(-1, 1)*delta_sigma); // new proposal for sigma
            }while(sigma_new<0.1);
            vector <double> H_new = compute_H( rnd, mu_new, sigma_new); // new proposal for H
            if (metro_par(rnd, H[0] , H_new[0] , beta)){
                mu = mu_new; // accept the new proposal
                sigma = sigma_new; // accept the new proposal
                H = H_new;
            } else {
                H = compute_H(rnd, mu, sigma); // update H anyway, even if not accepted
            }
        }
        coutp << T << "\t \t" << mu << "\t \t" << sigma << "\t \t " << H[0] << "\t \t" << H[1] << endl;

    }

    return 0;
}

vector<double> compute_H ( Random& rnd , double mu, double sigma ){
    // Simulation parameters
    double x = 0.0; // initial configuration (position)
    double delta = 0.5; // amplitude of Metropolis step

    int n_steps = 1000; // number of steps per block
    int n_blocks = 100; // number of blocks
    int tune_steps = 1000; // number of steps for delta tuning
    double tuned_acceptance = set_delta(rnd, delta, mu, sigma, x, tune_steps);
    
    double av{}; // when inside the for loop, contains average of averages of first (i+1) blocks
    double av2{}; // when inside the for loop, average of squared averages of first (i+1) blocks
    double sum_prog{}; // accumulates sum of averages of first (i+1) blocks
    double sum_prog2{}; // accumulates sum of squared averages

    vector<double> H;

    int block_eq = 0; // number of blocks discarded for equilibration

    for(int i=0; i <n_blocks; i++){ // loop over blocks
        double stima = 0;
        for(int j=0; j < n_steps; j++){ // loop over steps in a block
            // Metropolis step
            double y = x +  rnd.Rannyu(-1.0,1.0) * delta;
            if(metro_H(x,y,mu, sigma,rnd)){
                x = y;
            }
            if (i < block_eq)
                    continue;
            stima +=  Hpsi(x, mu, sigma);
        }
        sum_prog= sum_prog + stima/n_steps;
        sum_prog2 =sum_prog2 + stima/n_steps * stima/n_steps;
        av = (sum_prog) / (i+1);
        av2 = ( sum_prog2 ) / (i+1);
    }
    
  H.push_back(av);
  H.push_back(error(av, av2, n_blocks-1));

  return H;
}


double error( double AV, double AV2 , int n ){
    if(n == 0){
        return 0;
    }
    return sqrt( (AV2 - AV*AV)/n );
}

double psi(double x, double mu, double sigma){
    return exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
}

double Hpsi(double x , double mu, double sigma){
    return -0.5 *(  -1/(sigma*sigma) *  exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) -1/(sigma*sigma) *  exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) + 1/(pow(sigma,4)) * (x-mu)*(x-mu) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + 1/(pow(sigma,4)) * (x+mu)*(x+mu) * exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) )/psi(x, mu, sigma) + pow(x,4) - 5./2. * pow(x,2);
}

bool metro_H(double x , double y, double mu, double sigma, Random&  rnd){ // Metropolis algorithm
    bool decision = false;
    double acceptance = min( 1. , psi(y , mu , sigma)*psi(y , mu , sigma) / (psi(x, mu , sigma)*psi(x, mu,sigma)) ); 
    if(rnd.Rannyu() < acceptance ) decision = true; // Metropolis acceptance step
    return decision;
}

bool metro_par(Random& rnd, double H_old , double H_new , double beta ){
    bool decision = false;
    double acceptance = exp(-beta * ( H_new - H_old ) );
    if(rnd.Rannyu() < acceptance ) decision = true; // Metropolis acceptance step
    return decision;
}

double set_delta(Random& rnd, double& delta, double mu, double sigma, double x0, int steps) {
    double target_acceptance=0.5; // target value for acceptance
    double acceptance = 0;
    double x = x0; // starting point
    double tol = 0.01; // tolerance for acceptance close to target value
    int max_iter = 1000;
    int iter = 0;
    while (iter < max_iter) {
        int accepted = 0;
        for (int i = 0; i < steps; ++i) {
            double y = x + rnd.Rannyu(-1., 1.) * delta;
            if (metro_H(x, y, mu, sigma, rnd)) {
                x = y;
                accepted++;
            }
        }
        acceptance = static_cast<double>(accepted) / steps;
        if (fabs(acceptance - target_acceptance) < tol) {
            break;
        }
        if (acceptance < target_acceptance)
            delta *= 0.9; // reduce if acceptance too low
        else
            delta *= 1.1; // increase if too high

        iter++;
    }
    return acceptance;
}
