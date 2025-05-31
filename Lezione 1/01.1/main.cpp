#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "random.h" 

using namespace std;

// Funzione che calcola l'errore statistico
double error(vector<double> AV, vector<double> AV2, int n);

int main() {
// Parametri di simulazione
int M = 1000000;  // Numero totale di lanci casuali
int N = 100;      // Numero di blocchi
int L = M / N;    // Numero di lanci per blocco

Random rand;      // Oggetto per la generazione di numeri casuali
vector<double> av, av2; // Vettori per i valori medi e quadrati dei valori medi
vector<double> sum_prog(N), sum2_prog(N), err_prog(N); // Vettori per la somma progressiva e l'errore
ofstream out1, out2, out3; // File di output per i risultati

// Leggi il seme dalla sequenza di numeri primi
int seed[4];
int p1, p2;
ifstream Primes("Primes");
if (Primes.is_open()) {
    Primes >> p1 >> p2; // Legge due numeri primi
} else {
    cerr << "PROBLEM: Unable to open Primes" << endl; // Messaggio di errore se il file non si apre
}
Primes.close();

// Leggi il seme per la generazione di numeri casuali dal file seed.in
ifstream input("seed.in");
string property;
if (input.is_open()) {
    while (!input.eof()) {
        input >> property;
        if (property == "RANDOMSEED") {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rand.SetRandom(seed, p1, p2); // Inizializza il generatore di numeri casuali
        }
    }
    input.close();
} else {
    cerr << "PROBLEM: Unable to open seed.in" << endl; // Messaggio di errore se il file non si apre
}

// Calcola la media di L lanci casuali e la sua varianza (per la prima parte)
for (int i = 0; i < N; i++) {
    double sum1 = 0;
    for (int j = 0; j < L; j++) {
        sum1 += rand.Rannyu(); // Somma i numeri casuali generati
    }
    av.push_back(sum1 / L);   // Calcola la media per il blocco
    av2.push_back(pow(av[i], 2)); // Calcola il quadrato della media per il blocco
}

// Scrivi i risultati delle medie nei file di output
out1.open("media.dat");
out2.open("errore.dat");

for (int i = 0; i < N; i++) {
    // Somma progressiva per la media e il quadrato della media
    for (int j = 0; j <= i; j++) {
        sum_prog[i] += av[j];
        sum2_prog[i] += av2[j];
    }
    sum_prog[i] /= (i + 1); // Media progressiva
    out1 << sum_prog[i] << endl;
    sum2_prog[i] /= (i + 1); // Media progressiva del quadrato
    err_prog[i] = error(sum_prog, sum2_prog, i); // Calcola l'errore
    out2 << err_prog[i] << endl;
}

out1.close();
out2.close();

// Calcola la varianza e l'errore per la seconda parte
vector<double> av_, av2_; // Vettori per la seconda parte della simulazione
vector<double> sum_prog_(N), sum2_prog_(N), err_prog_(N);

for (int i = 0; i < N; i++) {
    double sum2 = 0;
    for (int j = 0; j < L; j++) {
        sum2 += pow(rand.Rannyu() - 0.5, 2); // Calcola la varianza di una variabile casuale uniformemente distribuita
    }
    av_.push_back(sum2 / L);   // Media della varianza per il blocco
    av2_.push_back(pow(av_[i], 2)); // Calcola il quadrato della media per la varianza
}

out1.open("var.dat");
out2.open("errorevar.dat");

for (int i = 0; i < N; i++) {
    // Somma progressiva per la media e il quadrato della media
    for (int j = 0; j <= i; j++) {
        sum_prog_[i] += av_[j];
        sum2_prog_[i] += av2_[j];
    }
    sum_prog_[i] /= (i + 1); // Media progressiva
    out1 << sum_prog_[i] << endl;
    sum2_prog_[i] /= (i + 1); // Media progressiva del quadrato
    err_prog_[i] = error(sum_prog_, sum2_prog_, i); // Calcola l'errore
    out2 << err_prog_[i] << endl;
}

out1.close();
out2.close();

// Calcola il chi-quadrato per una distribuzione uniforme
int m = 100; // Numero di intervalli
int n = 10000; // Numero di lanci
double expected = double(n) / m; // Valore atteso per ogni intervallo
vector<double> chi2(m); // Vettore per i valori di chi-quadrato
out3.open("chi2.dat");

for (int i = 0; i < m; i++) {
    vector<double> ni(m, 0); // Vettore per le frequenze osservate
    double chi2_ = 0;

    for (int j = 0; j < n; j++) {
        double x = rand.Rannyu(); // Genera un numero casuale
        int interval = floor(x * m); // Determina l'intervallo in cui cade il numero
        ni[interval]++; // Incrementa la frequenza osservata
    }

    // Calcola il chi-quadrato
    for (int j = 0; j < m; j++) {
        chi2_ += pow(ni[j] - expected, 2) / expected;
    }

    chi2[i] = chi2_; // Salva il valore di chi-quadrato
    out3 << chi2_ << endl;
}

out3.close();

// Salva il seme per poter riprodurre la simulazione
rand.SaveSeed();

return 0;
}

// Funzione che calcola l'errore statistico (deviazione standard)
double error(vector<double> AV, vector<double> AV2, int n) {
    if (n == 0) {
        return 0; // Se siamo al primo passo, l'errore è 0
    } else {
        return sqrt((AV2[n] - pow(AV[n], 2)) / n); // Calcolo dell'errore
    }
}

