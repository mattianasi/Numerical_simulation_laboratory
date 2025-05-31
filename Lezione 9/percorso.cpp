
#include <iostream>
#include <math.h>
#include "percorso.h"

using namespace std;

Route& Route::operator=(const Route& other) {
    if (this != &other) {
        this->_route = other._route;
        this->_length = other._length;
        this->_ndim = other._ndim;
        this->_point_rnd = other._point_rnd;
    }
    return *this;
}


void Route :: setstop(int stop , int city){
    _route(stop) = city;
    return;
}

int Route :: getstop(int stop){
    return _route(stop);
}


bool Route ::  check(){
    for(int i = 0; i < _ndim ; i++){
        for( int j = i+1; j < _ndim ; j++){
            if(_route(i) == _route(j)){
                return false;
            }
        }
    }
    return true;
}

int Route :: pbc(int city) const{
    if(city >= _ndim){
        return city - _ndim + 1;
    }
    return city;
}

void Route :: swap(int i, int j){
    int temp = _route(i);
    _route(i) = _route(j);
    _route(j) = temp;
    return;
}

void Route :: swap(){
    int i = (int) _point_rnd->Rannyu(1,_ndim); // generates a random number in [1,_ndim)
    int j = (int) _point_rnd->Rannyu(1,_ndim);
    int temp = _route(i);
    _route(i) = _route(j);
    _route(j) = temp;
    return;
}

double Route :: calculate_length( const arma::mat * distance_matrix ) const {
    double length = 0;
    for(int i=0 ; i < _ndim - 1 ; i++){
        length += distance_matrix->at(_route(i)-1, _route(pbc(i+1))-1);
    }
    length += distance_matrix->at(_route(_ndim-1)-1, _route(0)-1);
    //_length = length;
    return length;
}


void Route :: initialize(  arma::mat * distance_matrix , Random &rnd ){

    _point_rnd = &rnd;
    _route.resize(_ndim);

    for(int i=0 ; i< _ndim ; i ++){
      this->setstop(i,i+1);
    }

    int n_swaps = _ndim;
    for(int i=0 ; i < n_swaps ; i++){
        int j = (int) rnd.Rannyu(1,_ndim); // generates a random number in [1,_ndim)
        int k = (int) rnd.Rannyu(1,_ndim);
        this->swap(j,k);
    }

    std::ifstream infile("cities.dat");
    if (!infile.is_open()) {
        std::cerr << "Errore: impossibile aprire il file cities.dat" << std::endl;
        exit(1);
    }
    
    _length = this->calculate_length(distance_matrix);


   return;
}


/*void Route :: shift( ){
    int n = (int) _point_rnd->Rannyu(1,_ndim-1);
    arma::Col<int> temp(_ndim);
    for(int i=0 ; i< _ndim ; i++){
        temp(i) = _route(i);
    }
    for(int j = 1; j < _ndim ; j++){
        _route(pbc(j+n)) = temp(j);
    }
    if(check() == false){
        cout << "Error: Route not valid after shift" << endl;
    }
    return;
}*/

void Route::shift() {
    int m = (int) _point_rnd->Rannyu(1, _ndim - 2); // m < N-1
    int start = (int) _point_rnd->Rannyu(1, _ndim - m); // start from 1 to N - m - 1
    int n = (int) _point_rnd->Rannyu(1, _ndim - m - start + 1); // ensure shift doesn't go past end

    arma::Col<int> new_route(_ndim);
    int idx = 0;

    // Copia città fino al blocco da shiftare
    for(int i = 0; i < start; ++i)
        new_route(idx++) = _route(i);

    // Copia le città tra il blocco e la posizione di inserimento
    for(int i = start + m; i < start + m + n && i < _ndim; ++i)
        new_route(idx++) = _route(i);

    // Copia il blocco da shiftare
    for(int i = start; i < start + m; ++i)
        new_route(idx++) = _route(i);

    // Copia il resto
    for(int i = start + m + n; i < _ndim; ++i)
        new_route(idx++) = _route(i);

    _route = new_route;

    if (!check()) {
        cout << "Error: Route not valid after shift" << endl;
    }
}


/*void Route :: swap_block(){
    int m = (int) _point_rnd->Rannyu(1,_ndim/2 - 1 );
    int pos =(int) _point_rnd->Rannyu(1,_ndim - 1);

    for(int i=0 ; i < m ; i++){
        int j = pbc(i+pos);
        int k = pbc(i+m+pos);
        this->swap(j,k);
    }
    if(check() == false){
        cout << "Error: Route not valid after swap_block" << endl;
    }
}*/

void Route::swap_block() {
    int m = (int)_point_rnd->Rannyu(1, _ndim / 2);
    int pos1 = (int)_point_rnd->Rannyu(1, _ndim - 2 * m);
    int pos2 = (int)_point_rnd->Rannyu(pos1 + m, _ndim - m); // garantisce non sovrapposizione

    arma::Col<int> block1 = _route.subvec(pos1, pos1 + m - 1);
    arma::Col<int> block2 = _route.subvec(pos2, pos2 + m - 1);

    _route.subvec(pos1, pos1 + m - 1) = block2;
    _route.subvec(pos2, pos2 + m - 1) = block1;

    if(!check()) cout << "Error: Route not valid after swap_block" << endl;
}


void Route :: invert_block(){
    int start = pbc((int) _point_rnd->Rannyu(1, _ndim - 1 ));
    int end = pbc((int) _point_rnd->Rannyu(1,_ndim - 1));

    if (start > end) {
        std::swap(start, end);
    }
    arma::Col<int> temp(_ndim);
    for(int i=0 ; i< _ndim ; i++){
        temp(i) = _route(i);
    }

    while (start < end) {
        _route(start) = temp(end);
        _route(end) = temp(start);
        start++;
        end--;
    }
    if(check() == false){
        cout << "Error: Route not valid after invert_block" << endl;
    }
    return;

}


void Route :: mutate(){
    double i = _point_rnd->Rannyu();
    if(i < _pmut){
        double sel = _point_rnd->Rannyu();
        if(sel < 0.25){
            swap();
            return;
        } else if(sel < 0.5) {
            shift();
            return;
        } else if(sel < 0.75) {
            swap_block();
            return;
        } else {
            invert_block();
            return;
        }
    }
    return;
}