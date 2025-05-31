#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "random.h"

using namespace std;

double error(vector<double>, vector<double>, int);
double function(double);
double function1(double);

int main(){

Random rand;

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
            rand.SetRandom(seed,p1,p2);
        }
    }
    input.close();
} else cerr << "PROBLEM: Unable to open seed.in" << endl;

int M = 100000;
int N = 100;
int L = M/N;

vector<double> sum(N);
vector<double> sum2(N);
vector<double> err(N);
ofstream out1,out2;

for(int i=0; i<N; i++){
    double f = 0;
    for(int j=0; j<L; j++){
        double x = rand.Rannyu();
        f += function(x);
    }
    sum[i] = f/L;
    sum2[i] = sum[i]*sum[i];
}

out1.open("data.dat");
vector<double> averages(N, 0.0);
vector<double> averages2(N, 0.0);

for (int i = 0; i < N; i++) {
    for (int j = 0; j < i + 1; j++) {
        averages[i] += sum[j];
        averages2[i] += sum2[j];
    }
    averages[i] /= (i + 1);
    averages2[i] /= (i + 1);
    err[i] = error(averages, averages2, i);
    out1 << averages[i] << "," << err[i] << endl;
}
out1.close();

vector<double> sum_(N);
vector<double> sum2_(N);
vector<double> err_(N);

for(int i=0; i<N; i++){
    double f = 0;
    for(int j=0; j<L; j++){
        double x = rand.Distr();
        f += function1(x);
    }
    sum_[i] = f/L;
    sum2_[i] = sum_[i]*sum_[i];
}

out2.open("data1.dat");
vector<double> averages_(N, 0.0);
vector<double> averages2_(N, 0.0);

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

double error(vector<double> AV, vector<double> AV2, int n){
        if(n == 0) {
            return 0;}
        else{
            return sqrt((AV2[n]-pow(AV[n],2))/n);}
};

double function(double x){
    return M_PI/2*cos(M_PI*x/2);
}

double function1(double x){
    return M_PI/4*(cos(M_PI*x/2))/(1-x);
}