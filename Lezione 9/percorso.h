#ifndef __percorso__
#define __percorso__

#include <armadillo>
#include "random.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace arma;

// Class representing a single TSP route (individual in the population)
class Route {

private:
  int _ndim = 34;                        // Number of cities
  arma::Col<int> _route;                // Permutation of city indices (the path)
  double _length;                       // Total length of the route
  Random* _point_rnd;                   // Pointer to shared RNG
  double _pmut = 0.2;                   // Probability of applying a mutation

public:
  Route& operator=(const Route& other); // Assignment operator

  void initialize(arma::mat* distance_matrix, Random& rnd); // Generate an initial route
  void setstop(int stop, int city);                         // Set city index at given position
  int getstop(int stop);                                    // Get city index at given position
  bool check();                                             // Validate the route (e.g. all cities once)
  int pbc(int city) const;                                  // Handle periodic indexing
  double calculate_length(const arma::mat* distance_matrix) const; // Compute total route length
  void setlength(const arma::mat* distance_matrix);         // Update the cached length

  int getdim() { return _ndim; }                            // Return number of cities
  double getlength() { return _length; }                    // Return cached length
  arma::Col<int> getroute() { return _route; }              // Return the current route
  void setroute(arma::Col<int> route) { _route = route; }   // Set a new route manually

  // Mutation operators
  void swap(int i, int j);    // Swap cities at two positions
  void swap();                // Apply random pairwise swap
  void shift();               // Circularly shift a subsequence
  void swap_block();          // Swap two blocks of cities
  void invert_block();        // Invert the order of a segment
  void mutate();              // Apply a random mutation with probability _pmut
};

#endif // __percorso__
