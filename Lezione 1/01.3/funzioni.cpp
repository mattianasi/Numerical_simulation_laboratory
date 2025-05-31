#include "funzioni.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

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