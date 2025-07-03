#ifndef __percorso__
#define __percorso__

#include <armadillo>
#include "random.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace arma;

class Route {

private:
  int _ndim = 110; // Dimensionality of the system (number of cities)
  arma::Col<int> _route; // Vector storing the order of cities in the route
  double _length; // Total length of the current route
  Random* _point_rnd; // Pointer to random number generator
  double _pmut = 0.3; // Probability for performing mutation

public: 
  Route& operator=(const Route& other); 
  void initialize( arma::mat * distance_matrix , Random &rnd ); // Randomly initialize the route and compute its length
  void setstop(int stop , int city); // Set the city at a specific position in the route
  int getstop(int stop); // Get the city at a specific position in the route
  bool check(); // Check validity of the route (e.g., all cities visited once)
  int pbc(int city) const; // Apply periodic boundary conditions to city index
  double calculate_length( const arma::mat * distance_matrix ) const ; // Compute the total route length
  void setlength(const arma::mat * distance_matrix ); // Set internal _length variable from a distance matrix
  int getdim(){ return _ndim; }
  double getlength(){ return _length; }
  arma::Col<int> getroute(){ return _route; }
  void setroute(arma::Col<int> route){ _route = route; }

  // Mutation functions
  void swap(int i , int j); // Swap two cities at positions i and j
  void swap(); // Perform a random swap mutation
  void shift(); // Shift a subsequence of the route to a different position
  void swap_block(); // Swap two blocks (subsequences) of the route
  void invert_block(); // Invert the order of cities in a block
  void mutate(); // Apply a random mutation based on _pmut
};

#endif // __Route__
