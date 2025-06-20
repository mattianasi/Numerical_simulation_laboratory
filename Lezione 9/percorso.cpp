#include <iostream>
#include <math.h>
#include "percorso.h"

using namespace std;

// Assignment operator: copy all fields except random number generator
Route& Route::operator=(const Route& other) {
    if (this != &other) {
        this->_route = other._route;
        this->_length = other._length;
        this->_ndim = other._ndim;
        this->_point_rnd = other._point_rnd;
    }
    return *this;
}

// Set the city at a specific position in the route
void Route::setstop(int stop, int city) {
    _route(stop) = city;
}

// Get the city at a specific position
int Route::getstop(int stop) {
    return _route(stop);
}

// Check if route is valid (all cities visited once)
bool Route::check() {
    for (int i = 0; i < _ndim; i++) {
        for (int j = i + 1; j < _ndim; j++) {
            if (_route(i) == _route(j)) {
                return false;
            }
        }
    }
    return true;
}

// Periodic boundary condition helper
int Route::pbc(int city) const {
    if (city >= _ndim) {
        return city - _ndim + 1;
    }
    return city;
}

// Swap two cities in the route
void Route::swap(int i, int j) {
    std::swap(_route(i), _route(j));
}

// Randomly swap two cities (excluding the first)
void Route::swap() {
    int i = (int)_point_rnd->Rannyu(1, _ndim);
    int j = (int)_point_rnd->Rannyu(1, _ndim);
    std::swap(_route(i), _route(j));
}

// Compute total route length based on distance matrix
double Route::calculate_length(const arma::mat* distance_matrix) const {
    double length = 0;
    for (int i = 0; i < _ndim - 1; i++) {
        length += distance_matrix->at(_route(i) - 1, _route(pbc(i + 1)) - 1);
    }
    length += distance_matrix->at(_route(_ndim - 1) - 1, _route(0) - 1); // return to start
    return length;
}

// Initialize route with a random permutation and calculate its length
void Route::initialize(arma::mat* distance_matrix, Random& rnd) {
    _point_rnd = &rnd;
    _route.resize(_ndim);

    for (int i = 0; i < _ndim; i++) {
        _route(i) = i + 1;
    }

    // Shuffle the route with random swaps
    for (int i = 0; i < _ndim; i++) {
        int j = (int)rnd.Rannyu(1, _ndim);
        int k = (int)rnd.Rannyu(1, _ndim);
        this->swap(j, k);
    }

    // Optionally read city data if needed (currently just checks file exists)
    std::ifstream infile("cities.dat");
    if (!infile.is_open()) {
        std::cerr << "Errore: impossibile aprire il file cities.dat" << std::endl;
        exit(1);
    }

    _length = this->calculate_length(distance_matrix);
}

// Shift a block of cities forward by 'n' positions
void Route::shift() {
    int m = (int)_point_rnd->Rannyu(1, _ndim - 2);
    int start = (int)_point_rnd->Rannyu(1, _ndim - m);
    int n = (int)_point_rnd->Rannyu(1, _ndim - m - start + 1);

    arma::Col<int> new_route(_ndim);
    int idx = 0;

    // Copy cities before the shift block
    for (int i = 0; i < start; ++i)
        new_route(idx++) = _route(i);

    // Copy the gap between destination and source
    for (int i = start + m; i < start + m + n && i < _ndim; ++i)
        new_route(idx++) = _route(i);

    // Insert the block
    for (int i = start; i < start + m; ++i)
        new_route(idx++) = _route(i);

    // Copy remaining cities
    for (int i = start + m + n; i < _ndim; ++i)
        new_route(idx++) = _route(i);

    _route = new_route;

    if (!check()) {
        cout << "Error: Route not valid after shift" << endl;
    }
}

// Swap two non-overlapping blocks of equal size
void Route::swap_block() {
    int m = (int)_point_rnd->Rannyu(1, _ndim / 2);
    int pos1 = (int)_point_rnd->Rannyu(1, _ndim - 2 * m);
    int pos2 = (int)_point_rnd->Rannyu(pos1 + m, _ndim - m);

    arma::Col<int> block1 = _route.subvec(pos1, pos1 + m - 1);
    arma::Col<int> block2 = _route.subvec(pos2, pos2 + m - 1);

    _route.subvec(pos1, pos1 + m - 1) = block2;
    _route.subvec(pos2, pos2 + m - 1) = block1;

    if (!check()) cout << "Error: Route not valid after swap_block" << endl;
}

// Invert the order of cities between two random indices
void Route::invert_block() {
    int start = pbc((int)_point_rnd->Rannyu(1, _ndim - 1));
    int end = pbc((int)_point_rnd->Rannyu(1, _ndim - 1));
    if (start > end) std::swap(start, end);

    arma::Col<int> temp(_route);
    while (start < end) {
        _route(start) = temp(end);
        _route(end) = temp(start);
        ++start;
        --end;
    }

    if (!check()) {
        cout << "Error: Route not valid after invert_block" << endl;
    }
}

// Apply a random mutation based on probabilities
void Route::mutate() {
    double i = _point_rnd->Rannyu();
    if (i < _pmut) {
        double sel = _point_rnd->Rannyu();
        if (sel < 0.25) {
            swap();
        } else if (sel < 0.5) {
            shift();
        } else if (sel < 0.75) {
            swap_block();
        } else {
            invert_block();
        }
    }
}
