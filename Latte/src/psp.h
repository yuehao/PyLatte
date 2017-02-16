//
// Created by Yue Hao on 8/24/16.
//

#ifndef LATTE_PSP_H
#define LATTE_PSP_H
#include "tpsa.h"
#include <armadillo>
#include <vector>

class psp {
    //Power series 'particle'
private:
    int _dim;  //Number of phase space variable
    int _nvar;
    std::vector<CTPS> _pslist;
public:
    psp();
    CTPS &operator[](const int &i) { return _pslist[i]; }
    const CTPS &operator[](const int &i) const { return _pslist[i]; }
    void unitary();
    arma::mat getmap() const;
    void clear();
    std::vector<double> getcoor() const;
    
};


#endif //LATTE_PS6D_H
