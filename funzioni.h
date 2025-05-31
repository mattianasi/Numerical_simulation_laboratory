#ifndef __funzioni__
#define __funzioni__

#include "random.h"

void LoadPrimes(int &p1, int &p2);
void LoadSeed(int seed[4], Random &rand, int p1, int p2);
double error(const vector<double>& averages, const vector<double>& averages2, int i);
double CalculateDistance(const vector<int>& pos);




#endif