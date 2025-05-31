#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "funzioni.h"

using namespace std;

double error(vector<double>, vector<double>, int);

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


int M = 100000, N = 100, L = M / N; 
double S0 = 100, T = 1, K=100, r=0.1, sigma=0.25;

// Vettori per memorizzare i risultati
vector<double> av(N), av2(N), err(N);

// Loop su ciascun esperimento
for(int i = 0; i < N; i++){
    double C=0;
    for(int j = 0; j < L; j++){
        double z = rand.Gauss(0,1);
        double S = S0*exp((r-sigma*sigma/2)*T+sigma*z*sqrt(T));
        C += exp(-r*T)*fmax(0,S-K);
    }
    av[i] = C/L;  // Calcola la stima di pi per ogni esperimento
    av2[i] = pow(av[i], 2);  // Per il calcolo dell'errore
}
// Scrittura dei risultati su file
ofstream out("data1.dat");
vector<double> sum_prog(N, 0.0), sum_prog2(N, 0.0);

// Calcolo delle medie e dell'errore per ogni esperimento
for (int i = 0; i < N; i++) {
    for (int j = 0; j <= i; j++) {
        sum_prog[i] += av[j];  // Somma delle stime di pi
        sum_prog2[i] += av2[j];  // Somma dei quadrati delle stime di pi
    }
    sum_prog[i] /= (i + 1);  // Media delle stime fino all'indice i
    sum_prog2[i] /= (i + 1);  // Media dei quadrati delle stime fino all'indice i
    err[i] = error(sum_prog, sum_prog2, i); // Calcola l'errore 
    out << sum_prog[i] << "\t" << err[i] << endl;  // Scrive i risultati nel file
}
out.close();

//

for(int i = 0; i < N; i++){
    double P=0;
    for(int j = 0; j < L; j++){
        double z = rand.Gauss(0,1);
        double S = S0*exp((r-sigma*sigma/2)*T+sigma*z*sqrt(T));
        P += exp(-r*T)*fmax(0,K-S);
    }
    av[i] = P/L;  // Calcola la stima di pi per ogni esperimento
    av2[i] = pow(av[i], 2);  // Per il calcolo dell'errore
}
// Scrittura dei risultati su file
out.open("data2.dat");

// Calcolo delle medie e dell'errore per ogni esperimento
for (int i = 0; i < N; i++) {
    sum_prog[i] = 0;
    sum_prog2[i] = 0;
    err[i] = 0;
    for (int j = 0; j <= i; j++) {
        sum_prog[i] += av[j];  // Somma delle stime di pi
        sum_prog2[i] += av2[j];  // Somma dei quadrati delle stime di pi
    }
    sum_prog[i] /= (i + 1);  // Media delle stime fino all'indice i
    sum_prog2[i] /= (i + 1);  // Media dei quadrati delle stime fino all'indice i
    err[i] = error(sum_prog, sum_prog2, i); // Calcola l'errore 
    out << sum_prog[i] << "\t" << err[i] << endl;  // Scrive i risultati nel file
}
out.close();

//

for(int i = 0; i < N; i++){
    double C=0;
    for(int j = 0; j < L; j++){
        double S1 = S0;
        double S2;
        for(int k=0; k<100; k++){
            double z = rand.Gauss(0,1);
            S2 = S1*exp((r-sigma*sigma/2)*1.0/100+sigma*z*sqrt(1.0/100));
            S1 = S2;
        }
        C += exp(-r*T)*fmax(0,S2-K);
    }
    av[i] = C/L;  // Calcola la stima di pi per ogni esperimento
    av2[i] = pow(av[i], 2);  // Per il calcolo dell'errore
}
// Scrittura dei risultati su file
out.open("data3.dat");

// Calcolo delle medie e dell'errore per ogni esperimento
for (int i = 0; i < N; i++) {
    sum_prog[i] = 0;
    sum_prog2[i] = 0;
    err[i] = 0;
    for (int j = 0; j <= i; j++) {
        sum_prog[i] += av[j];  // Somma delle stime di pi
        sum_prog2[i] += av2[j];  // Somma dei quadrati delle stime di pi
    }
    sum_prog[i] /= (i + 1);  // Media delle stime fino all'indice i
    sum_prog2[i] /= (i + 1);  // Media dei quadrati delle stime fino all'indice i
    err[i] = error(sum_prog, sum_prog2, i); // Calcola l'errore 
    out << sum_prog[i] << "\t" << err[i] << endl;  // Scrive i risultati nel file
}
out.close();

//

for(int i = 0; i < N; i++){
    double P=0;
    for(int j = 0; j < L; j++){
        double S1 = S0;
        double S2;
        for(int k=0; k<100; k++){
            double z = rand.Gauss(0,1);
            S2 = S1*exp((r-sigma*sigma/2)*1.0/100+sigma*z*sqrt(1.0/100));
            S1 = S2;
        }
        P += exp(-r*T)*fmax(0,K-S2);
    }
    av[i] = P/L;  // Calcola la stima di pi per ogni esperimento
    av2[i] = pow(av[i], 2);  // Per il calcolo dell'errore
}
// Scrittura dei risultati su file
out.open("data4.dat");

// Calcolo delle medie e dell'errore per ogni esperimento
for (int i = 0; i < N; i++) {
    sum_prog[i] = 0;
    sum_prog2[i] = 0;
    err[i] = 0;
    for (int j = 0; j <= i; j++) {
        sum_prog[i] += av[j];  // Somma delle stime di pi
        sum_prog2[i] += av2[j];  // Somma dei quadrati delle stime di pi
    }
    sum_prog[i] /= (i + 1);  // Media delle stime fino all'indice i
    sum_prog2[i] /= (i + 1);  // Media dei quadrati delle stime fino all'indice i
    err[i] = error(sum_prog, sum_prog2, i); // Calcola l'errore 
    out << sum_prog[i] << "\t" << err[i] << endl;  // Scrive i risultati nel file
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