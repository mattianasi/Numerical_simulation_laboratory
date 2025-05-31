#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"

using namespace std;

int main(){

Random rand;
int N = 10000;
vector<double> S1(N);
vector<double> S2(N);
vector<double> S10(N);
vector<double> S100(N);
ofstream out1, out2, out3;

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

out1.open("lin.dat");
for(int i=0; i<N; i++){
    S1[i] = rand.Rannyu();
    S2[i] = (rand.Rannyu()+rand.Rannyu())/2;
    for(int j=0; j<10; j++){
        S10[i] += rand.Rannyu();
    }
    for(int j=0; j<100; j++){
        S100[i] += rand.Rannyu();
    }
    out1 << S1[i] << "," << S2[i] << "," << S10[i]/10 << "," << S100[i]/100 << endl;
}
out1.close();

out2.open("exp.dat");
for(int i=0; i<N; i++){
    S1[i] = rand.Exp(1);
    S2[i] = (rand.Exp(1)+rand.Exp(1))/2;
    for(int j=0; j<10; j++){
        S10[i] += rand.Exp(1);
    }
    for(int j=0; j<100; j++){
        S100[i] += rand.Exp(1);
    }
    out2 << S1[i] << "," << S2[i] << "," << S10[i]/10 << "," << S100[i]/100 << endl;
}
out2.close();

out3.open("lorentz.dat");
for(int i=0; i<N; i++){
    S1[i] = rand.Lorentz(0,1);
    S2[i] = (rand.Lorentz(0,1)+rand.Lorentz(0,1))/2;
    for(int j=0; j<10; j++){
        S10[i] += rand.Lorentz(0,1);
    }
    for(int j=0; j<100; j++){
        S100[i] += rand.Lorentz(0,1);
    }
    out3 << S1[i] << "," << S2[i] << "," << S10[i]/10 << "," << S100[i]/100 << endl;
}
out3.close();

return 0;
}