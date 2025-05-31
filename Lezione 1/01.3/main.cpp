#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "funzioni.h"

using namespace std;

double error(vector<double> AV, vector<double> AV2, int n);
bool Intersect(double l, double d, double ang, Random& rnd);

int main(){
Random rand;
int seed[4];
int p1, p2;

// Lettura dei numeri primi per inizializzare il generatore di numeri casuali
ifstream Primes("Primes");
if (Primes.is_open()){
Primes >> p1 >> p2;
} else cerr << "PROBLEM: Unable to open Primes" << endl;
Primes.close();

// Lettura del file di seed iniziale
ifstream input("seed.in");
string property;
if (input.is_open()){
while (!input.eof()){
    input >> property;
    if (property == "RANDOMSEED"){
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rand.SetRandom(seed, p1, p2);
    }
}
input.close();
} else cerr << "PROBLEM: Unable to open seed.in" << endl;

// ~~~~ Calcolo di pi greco con il metodo di Buffon ~~~~ //

int M = 1000000, N = 100, L = M / N, l = 1;  // M = numero di lanci, N = numero di esperimenti
double d = 3;  // Lunghezza della distanza tra le linee

// Vettori per memorizzare i risultati
vector<double> sum_pi(N), sum_pi2(N), err_pi(N);
double theta;

// Loop su ciascun esperimento
for(int i = 0; i < N; i++){
    int Nhit = 0;  // Numero di intersezioni per esperimento
    for(int j = 0; j < L; j++){
        double x = rand.Rannyu(-1,1);
        double y = rand.Rannyu(-1,1);
        if(y>=0){
            theta = acos(x/(sqrt(x*x+y*y)));
        }else{
            theta = 2*M_PI - acos(x/(sqrt(x*x+y*y)));
        }
        if(Intersect(l, d, theta/2, rand)) 
        Nhit++;  // Conta il numero di intersezioni
    }
    sum_pi[i] = (2.0 * l * L) / (d * Nhit);  // Calcola la stima di pi per ogni esperimento
    sum_pi2[i] = pow(sum_pi[i], 2);  // Per il calcolo dell'errore
}
// Scrittura dei risultati su file
ofstream out("pi.dat");
vector<double> sum_prog(N, 0.0), sum_prog2(N, 0.0);

// Calcolo delle medie e dell'errore per ogni esperimento
for (int i = 0; i < N; i++) {
    for (int j = 0; j <= i; j++) {
        sum_prog[i] += sum_pi[j];  // Somma delle stime di pi
        sum_prog2[i] += sum_pi2[j];  // Somma dei quadrati delle stime di pi
    }
    sum_prog[i] /= (i + 1);  // Media delle stime fino all'indice i
    sum_prog2[i] /= (i + 1);  // Media dei quadrati delle stime fino all'indice i
    err_pi[i] = error(sum_prog, sum_prog2, i); // Calcola l'errore 
    out << sum_prog[i] << "\t" << err_pi[i] << endl;  // Scrive i risultati nel file
}
out.close();
cout << "end program" << endl;
return 0;
}

double error(vector<double> AV, vector<double> AV2, int n){
        if(n == 0) {
            return 0;}
        else{
            return sqrt((AV2[n]-pow(AV[n],2))/n);}
};

bool Intersect(double l, double d, double ang, Random& rnd){
    double y = (l/2) * sin(ang); // Calcola la coordinata y del punto di intersezione
    double x = rnd.Rannyu() * d; // Genera un punto casuale lungo la distanza d
    return y >= min(x, d - x); // Se la coordinata y >= alla distanza tra x e d ==> intersezione avviene
};