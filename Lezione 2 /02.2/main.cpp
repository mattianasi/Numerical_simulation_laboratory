#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>
#include "random.h"

using namespace std;

double error(vector<double>, vector<double>, int);


int main(){

Random rand;
int N = 100;
int M = 100000;
int L = M/N;
double x;
double step = 100;

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

vector<double> sum(N);
vector<double> sum2(N);
vector<double> err(N);
ofstream out1;
double interval[6] = {1.0/6, 2.0/6, 3.0/6, 4.0/6, 5.0/6, 1.0};

for(int i=0; i<N; i++){
    vector<double> v(N, 0.0);
    for(int j=0; j<L; j++){
        vector<double> pos(6, 0.0);
        for(int k = 0; k<step; k++){
            x = rand.Rannyu();
            for(int m=0; m<6; m++){
                if(x<=interval[m]){
                pos[m]+=1;
                break;
                }
            }
            v[k] += pow((pos[0]-pos[1]),2)+pow((pos[2]-pos[3]),2)+pow((pos[4]-pos[5]),2);
        }
    }
    for(int k = 0; k<step; k++){
        v[k]/=L;
        v[k] = sqrt(v[k]);
    }
    for(int k=0; k<step; k++){
        sum[k] += v[k];
        sum2[k] += v[k]*v[k];
    }
}

out1.open("data.dat");


for (int i = 0; i < step; i++) { //anche sui blocchi? in caso modificare
    sum[i]/=(step);
    sum2[i]/=(step);
    err[i] = sqrt((sum2[i]-pow(sum[i],2))/step); 
    out1 << sum[i] << "," << err[i] << endl;
}
out1.close();

vector<double> sum_cont(N);
vector<double> sum_cont2(N);
vector<double> error_cont(N);
ofstream out2;

for(int i=0; i<N; i++){
    vector<double> w(N, 0.0);
    for(int j=0; j<L; j++){
        vector<double> r(3,0.0);
        for(int k = 0; k<step; k++){
            double a = rand.Rannyu();
            double phi = rand.Rannyu(0,2*M_PI);
            double theta = acos(1-2*a);
            r[0] += sin(theta)*cos(phi);
            r[1] += sin(theta)*sin(phi);
            r[2] += cos(theta);
            w[k] += r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
        }
    }
    for(int k = 0; k<step; k++){
        w[k]/=L;
        w[k] = sqrt(w[k]);
    }
    for(int k=0; k<step; k++){
        sum_cont[k] += w[k];
        sum_cont2[k] += w[k]*w[k];
    }
}

out2.open("datacont.dat");


for (int i = 0; i < N; i++) {
    sum_cont[i]/=(step);
    sum_cont2[i]/=(step);
    error_cont[i] = sqrt((sum_cont2[i]-pow(sum_cont[i],2))/step);
    out2 << sum_cont[i] << "," << error_cont[i] << endl;
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