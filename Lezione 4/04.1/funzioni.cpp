#include "funzioni.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

void SeedSet(Random rand){
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
}

double error(vector<double> AV, vector<double> AV2, int n){
        if(n == 0) {
            return 0;}
        else{
            return sqrt((AV2[n]-pow(AV[n],2))/n);}
};

//01.3
bool Intersect(double l, double d, double ang, Random& rnd){
    double y = (l/2) * sin(ang); // Calcola la coordinata y del punto di intersezione
    double x = rnd.Rannyu() * d; // Genera un punto casuale lungo la distanza d
    return y >= min(x, d - x); // Se la coordinata y >= alla distanza tra x e d ==> intersezione avviene
};