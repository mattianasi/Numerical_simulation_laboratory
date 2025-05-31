#include "funzioni.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

void SeedSet(int i, Random rand){
    int seed[i];
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
}